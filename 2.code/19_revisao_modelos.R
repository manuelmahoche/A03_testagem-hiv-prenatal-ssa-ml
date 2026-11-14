# ============================================================================
# P06 - Testagem de HIV no pré-natal na África Subsaariana
# Etapa 19 - Respostas à revisão interna: modelos, incerteza, equidade e utilidade
# Autor: Manuel Mahoche
# Entrada: ds_clean.parquet, particoes.parquet, 06_modelos_finais.rds (hiperparâmetros),
#          pred_teste_interno.parquet
# Saídas (3.output/tabs/19_*.csv; figura 3.output/figs/19_fig_dca.png):
#   19_auc_ic.csv                 AUC com IC95% (bootstrap por PSU): interno, LOCO, Estrato B, 2015-17
#   19_agregado_vs_pais.csv       AUC agregada vs média por país (Estrato B e 2015-17)
#   19_brier_escalonado.csv       Brier, Brier nulo e Brier skill score (interno, B bruto)
#   19_recalibracao_cv.csv        recalibração AJUSTADA e AVALIADA em amostras distintas (CV por PSU)
#   19_recalibracao_cv_pais.csv   idem por país
#   19_sensibilidade_preditores.csv  modelo sem conhecimento PTV/estigma e sem bloco de informação
#   19_calibracao_subgrupo.csv    calibração por subgrupo (intercepto, inclinação, O/E, Brier)
#   19_limiares_subgrupo.csv      sensibilidade/especificidade/VPP/VPN por limiar e subgrupo
#   19_dca.csv                    benefício líquido (decision curve analysis)
#   19_auc_regiao_xgb.csv         AUC por região SSA, gradient boosting e random forest
#   19_tamanho_amostra.csv        critérios de Riley et al. (2019) para o tamanho de amostra
#   19_testagem_sem_cpn.csv       cobertura entre mulheres sem pré-natal
# Convenções: o modelo prevê p = P(testada). Para focalização usa-se risco r = 1 - p
# (probabilidade de NÃO ser testada); o "evento" de interesse é não ter sido testada.
# Métricas de predição não ponderadas (como nas etapas 06-07).
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(here); library(xgboost); library(ggplot2)
})
source(here("2.code", "R", "config.R")); source(here("2.code", "R", "funcoes.R"))
source(here("2.code", "R", "tema_figuras.R"))
set.seed(semente)
NB <- 1000   # réplicas de bootstrap

pred_all <- setdiff(vars_preditores, "regiao_id")
dc <- ler_parquet(file.path(path_process, "ds_clean.parquet")) |>
  select(id_mulher, pais, nome_pais, regiao_ssa, amostra, estrato, psu_id, riqueza, escolaridade, residencia,
         idade_grp, anc_qualquer, w_dhs, w_igual, all_of(pred_all), y = hiv_anc_resultado)
pa <- ler_parquet(file.path(path_process, "particoes.parquet")) |> select(id_mulher, particao, cv_fold)
dc <- left_join(dc, pa, by = "id_mulher")
for (v in pred_all) { x <- as.factor(dc[[v]]); if (anyNA(x)) { x <- addNA(x); levels(x)[is.na(levels(x))] <- "(sem dado)" }; dc[[v]] <- x }
y <- dc$y
idx_tr   <- which(dc$amostra == "principal" & dc$particao == "treino")
idx_te   <- which(dc$amostra == "principal" & dc$particao == "teste")
idx_A    <- which(dc$amostra == "principal" & dc$estrato == "A")
idx_B    <- which(dc$amostra == "principal" & dc$estrato == "B")
idx_S    <- which(dc$amostra == "sensibilidade")
hp <- readRDS(file.path("3.output", "mod", "06_modelos_finais.rds"))$hp_xgb

# ---------------------------------------------------------------------------
# Utilitários
# ---------------------------------------------------------------------------
auc_fast <- function(yy, p) {                      # AUC por postos (Mann-Whitney)
  n1 <- sum(yy == 1); n0 <- length(yy) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(p)[yy == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
clamp <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
cal_int <- function(yy, p) unname(coef(glm(yy ~ offset(qlogis(clamp(p))), family = binomial))[1])
cal_slope <- function(yy, p) unname(coef(glm(yy ~ qlogis(clamp(p)), family = binomial))[2])
brier <- function(yy, p) mean((p - yy)^2)
brier_nulo <- function(yy) mean(yy) * (1 - mean(yy))
bss <- function(yy, p) 1 - brier(yy, p) / brier_nulo(yy)
auprc <- function(yy, p) { o <- order(p, decreasing = TRUE); yy <- yy[o]; tp <- cumsum(yy); fp <- cumsum(1 - yy)
  sum((tp / (tp + fp)) * diff(c(0, tp / sum(yy)))) }

build_X <- function(pred) {
  X <- model.matrix(~ . - 1, data = dc[pred])
  X[, colSums(X[idx_A, , drop = FALSE]) > 0 | colSums(X[idx_B, , drop = FALSE]) > 0, drop = FALSE]
}
fit_xgb <- function(X, rows) {
  m <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = hp$max_depth,
                               eta = hp$eta, subsample = 0.8, colsample_bytree = 0.8, nthread = 0),
                 data = xgb.DMatrix(X[rows, , drop = FALSE], label = y[rows]), nrounds = hp$nrounds, verbose = 0)
  function(nr) as.numeric(predict(m, X[nr, , drop = FALSE]))
}

# Gera as predições de uma variante de preditores (interno, LOCO, B, 2015-17)
run_variante <- function(pred) {
  X <- build_X(pred)
  f_in  <- fit_xgb(X, idx_tr)
  interno <- dc[idx_te, c("pais", "psu_id")] |> mutate(y = y[idx_te], p = f_in(idx_te))
  loco <- bind_rows(lapply(sort(unique(dc$pais[idx_A])), function(pc) {
    tr <- idx_A[dc$pais[idx_A] != pc]; te <- idx_A[dc$pais[idx_A] == pc]
    dc[te, c("pais", "psu_id")] |> mutate(y = y[te], p = fit_xgb(X, tr)(te))
  }))
  f_A <- fit_xgb(X, idx_A)
  list(interno = interno, loco = loco,
       B = dc[idx_B, c("pais", "psu_id")] |> mutate(y = y[idx_B], p = f_A(idx_B)),
       sens = dc[idx_S, c("pais", "psu_id")] |> mutate(y = y[idx_S], p = f_A(idx_S)),
       ncol = ncol(X))
}

# Bootstrap por conglomerado (PSU) dentro de cada país
boot_conj <- function(df, B = NB) {
  paises <- sort(unique(df$pais))
  psu_rows <- lapply(paises, function(pc) { d <- which(df$pais == pc); split(d, df$psu_id[d]) })
  names(psu_rows) <- paises
  res_c <- matrix(NA_real_, B, length(paises), dimnames = list(NULL, paises)); res_pool <- numeric(B)
  for (b in seq_len(B)) {
    all_idx <- vector("list", length(paises))
    for (j in seq_along(paises)) {
      pr <- psu_rows[[j]]; s <- sample.int(length(pr), replace = TRUE)
      idx <- unlist(pr[s], use.names = FALSE); all_idx[[j]] <- idx
      res_c[b, j] <- auc_fast(df$y[idx], df$p[idx])
    }
    ia <- unlist(all_idx, use.names = FALSE); res_pool[b] <- auc_fast(df$y[ia], df$p[ia])
  }
  list(por_pais = res_c, media = rowMeans(res_c, na.rm = TRUE), agregado = res_pool)
}
q <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
resumo_auc <- function(df, escopo, modelo) {
  bt <- boot_conj(df)
  pc <- sort(unique(df$pais))
  est_c <- sapply(pc, function(k) auc_fast(df$y[df$pais == k], df$p[df$pais == k]))
  linhas <- tibble(modelo = modelo, escopo = escopo, nivel = pc, n = as.integer(table(df$pais)[pc]),
                   auc = est_c, lo = apply(bt$por_pais, 2, function(z) q(z)[1]), hi = apply(bt$por_pais, 2, function(z) q(z)[2]))
  est_pool <- auc_fast(df$y, df$p); est_med <- mean(est_c)
  bind_rows(linhas,
    tibble(modelo = modelo, escopo = escopo, nivel = "Média por país", n = nrow(df), auc = est_med, lo = q(bt$media)[1], hi = q(bt$media)[2]),
    tibble(modelo = modelo, escopo = escopo, nivel = "Agregado (países empilhados)", n = nrow(df), auc = est_pool, lo = q(bt$agregado)[1], hi = q(bt$agregado)[2]),
    tibble(modelo = modelo, escopo = escopo, nivel = "Média por país ponderada por n", n = nrow(df),
           auc = weighted.mean(est_c, linhas$n), lo = NA_real_, hi = NA_real_))
}

# ---------------------------------------------------------------------------
# 1) Variantes de preditores: modelo completo, sem conhecimento PTV + estigma, sem bloco de informação
# ---------------------------------------------------------------------------
variantes <- list(
  completo        = pred_all,
  sem_ptv_estigma = setdiff(pred_all, c("conhece_ptv", "estigma")),
  sem_informacao  = setdiff(pred_all, vars_bloco_info))
cat("Ajustando variantes (interno + LOCO + B + 2015-17)...\n")
res <- list()
for (nm in names(variantes)) { cat("  ", nm, "\n"); res[[nm]] <- run_variante(variantes[[nm]]) }

# Verificação de reprodução do modelo completo (deve bater com 06)
chk <- c(interno = auc_fast(res$completo$interno$y, res$completo$interno$p),
         loco = mean(sapply(split(res$completo$loco, res$completo$loco$pais), function(d) auc_fast(d$y, d$p))))
cat(sprintf("  Reprodução (completo): AUC interno %.3f (06: 0,792); LOCO média %.3f (06: 0,768)\n", chk[1], chk[2]))

# ---------------------------------------------------------------------------
# 2) AUC com IC95% (bootstrap por PSU) e agregado vs média por país
# ---------------------------------------------------------------------------
cat("Bootstrap por PSU (", NB, " réplicas)...\n", sep = "")
auc_tab <- bind_rows(lapply(names(res), function(nm) {
  r <- res[[nm]]
  bt_in <- boot_conj(r$interno)
  bind_rows(
    tibble(modelo = nm, escopo = "Interno (teste, Estrato A)", nivel = "Agregado", n = nrow(r$interno),
           auc = auc_fast(r$interno$y, r$interno$p), lo = q(bt_in$agregado)[1], hi = q(bt_in$agregado)[2]),
    resumo_auc(r$loco, "LOCO (16 países)", nm),
    resumo_auc(r$B, "Externo (Estrato B)", nm),
    resumo_auc(r$sens, "Janela 2015-17", nm))
})) |> mutate(across(c(auc, lo, hi), ~ round(.x, 3)))
write_csv(auc_tab, file.path(path_tabs, "19_auc_ic.csv"))

agr <- auc_tab |> filter(modelo == "completo", nivel %in% c("Média por país", "Agregado (países empilhados)", "Média por país ponderada por n")) |>
  select(escopo, nivel, auc, lo, hi)
# 2015-17 sem o Benim (única com cobertura baixa)
s_alt <- res$completo$sens |> filter(pais != "BJ")
agr <- bind_rows(agr, tibble(escopo = "Janela 2015-17 sem Benim", nivel = "Média por país",
                             auc = round(mean(sapply(split(s_alt, s_alt$pais), function(d) auc_fast(d$y, d$p))), 3), lo = NA, hi = NA),
                 tibble(escopo = "Janela 2015-17 sem Benim", nivel = "Agregado (países empilhados)",
                        auc = round(auc_fast(s_alt$y, s_alt$p), 3), lo = NA, hi = NA))
write_csv(agr, file.path(path_tabs, "19_agregado_vs_pais.csv"))

# ---------------------------------------------------------------------------
# 3) Brier escalonado (interno e Estrato B bruto)
# ---------------------------------------------------------------------------
b_in <- res$completo$interno; b_B <- res$completo$B
bs_tab <- bind_rows(
  tibble(amostra = "Interno (teste, Estrato A)", prev_tested = mean(b_in$y), brier = brier(b_in$y, b_in$p),
         brier_nulo = brier_nulo(b_in$y), bss = bss(b_in$y, b_in$p)),
  tibble(amostra = "Estrato B, modelo bruto", prev_tested = mean(b_B$y), brier = brier(b_B$y, b_B$p),
         brier_nulo = brier_nulo(b_B$y), bss = bss(b_B$y, b_B$p)))
write_csv(bs_tab |> mutate(across(where(is.numeric), ~ round(.x, 4))), file.path(path_tabs, "19_brier_escalonado.csv"))

# ---------------------------------------------------------------------------
# 4) Recalibração ajustada e avaliada em amostras distintas (CV por PSU, 5 folds, 20 repetições)
# ---------------------------------------------------------------------------
recal_cv <- function(df, tipo = c("platt", "intercepto"), por_pais = FALSE, K = 5, R = 20) {
  tipo <- match.arg(tipo)
  lp <- qlogis(clamp(df$p))
  out <- vector("list", R)
  for (r in seq_len(R)) {
    pc <- numeric(nrow(df)); grupos <- if (por_pais) split(seq_len(nrow(df)), df$pais) else list(seq_len(nrow(df)))
    for (g in grupos) {
      psus <- unique(paste(df$pais[g], df$psu_id[g])); fold <- setNames(sample(rep_len(seq_len(K), length(psus))), psus)
      f <- fold[paste(df$pais[g], df$psu_id[g])]
      for (k in seq_len(K)) {
        tr <- g[f != k]; te <- g[f == k]
        m <- if (tipo == "platt") glm(df$y[tr] ~ lp[tr], family = binomial)
             else glm(df$y[tr] ~ offset(lp[tr]), family = binomial)
        pc[te] <- if (tipo == "platt") plogis(coef(m)[1] + coef(m)[2] * lp[te]) else plogis(coef(m)[1] + lp[te])
      }
    }
    out[[r]] <- pc
  }
  out
}
resumo_recal <- function(df, preds, rotulo) {
  met <- sapply(preds, function(p) c(brier = brier(df$y, p), bss = bss(df$y, p), auc = auc_fast(df$y, p),
                                     cal_int = cal_int(df$y, p), cal_slope = cal_slope(df$y, p)))
  # BSS contra o nulo específico do país (cobertura conhecida de cada país)
  nulo_pais <- mean(unlist(lapply(split(df, df$pais), function(d) rep(brier_nulo(d$y), nrow(d)))))
  tibble(versao = rotulo, brier = mean(met["brier", ]), brier_nulo_agregado = brier_nulo(df$y),
         bss_agregado = mean(met["bss", ]), brier_nulo_pais = nulo_pais,
         bss_nulo_pais = mean(1 - met["brier", ] / nulo_pais), auc = mean(met["auc", ]),
         cal_int = mean(met["cal_int", ]), cal_slope = mean(met["cal_slope", ]))
}
cat("Recalibração com validação cruzada por PSU (Estrato B)...\n")
dfB <- res$completo$B
rc <- bind_rows(
  resumo_recal(dfB, list(dfB$p), "Bruto (modelo do Estrato A)"),
  { lpB <- qlogis(clamp(dfB$p)); m <- glm(dfB$y ~ lpB, family = binomial)
    resumo_recal(dfB, list(as.numeric(fitted(m))), "Platt ajustado e avaliado na MESMA amostra (como na v1)") },
  resumo_recal(dfB, recal_cv(dfB, "intercepto"), "Só intercepto, agregado, CV por PSU"),
  resumo_recal(dfB, recal_cv(dfB, "platt"), "Platt (intercepto + inclinação), agregado, CV por PSU"),
  resumo_recal(dfB, recal_cv(dfB, "intercepto", por_pais = TRUE), "Só intercepto, por país, CV por PSU"),
  resumo_recal(dfB, recal_cv(dfB, "platt", por_pais = TRUE), "Platt, por país, CV por PSU"))
write_csv(rc |> mutate(across(where(is.numeric), ~ round(.x, 4))), file.path(path_tabs, "19_recalibracao_cv.csv"))

rc_pais <- bind_rows(lapply(sort(unique(dfB$pais)), function(pc) {
  d <- dfB[dfB$pais == pc, ]; pl <- recal_cv(d, "platt", por_pais = TRUE)
  resumo_recal(d, pl, "Platt por país, CV por PSU") |> mutate(pais = pc, n = nrow(d), cobertura_pct = round(100 * mean(d$y), 1),
                                                              brier_bruto = brier(d$y, d$p), bss_bruto = bss(d$y, d$p))
})) |> select(pais, n, cobertura_pct, brier_bruto, bss_bruto, brier, brier_nulo_agregado, bss_agregado, auc, cal_int, cal_slope)
write_csv(rc_pais |> mutate(across(where(is.numeric), ~ round(.x, 4))), file.path(path_tabs, "19_recalibracao_cv_pais.csv"))

# ---------------------------------------------------------------------------
# 5) Sensibilidade aos preditores medidos depois do desfecho
# ---------------------------------------------------------------------------
sens_tab <- bind_rows(lapply(names(res), function(nm) {
  r <- res[[nm]]
  loco_pc <- sapply(split(r$loco, r$loco$pais), function(d) auc_fast(d$y, d$p))
  B_pc <- sapply(split(r$B, r$B$pais), function(d) auc_fast(d$y, d$p))
  S_pc <- sapply(split(r$sens, r$sens$pais), function(d) auc_fast(d$y, d$p))
  tibble(variante = nm, n_colunas = r$ncol, auc_interno = auc_fast(r$interno$y, r$interno$p),
         auprc_interno = auprc(r$interno$y, r$interno$p), brier_interno = brier(r$interno$y, r$interno$p),
         bss_interno = bss(r$interno$y, r$interno$p), cal_int = cal_int(r$interno$y, r$interno$p),
         cal_slope = cal_slope(r$interno$y, r$interno$p), auc_loco_media = mean(loco_pc),
         auc_loco_min = min(loco_pc), auc_loco_max = max(loco_pc), auc_B_media = mean(B_pc),
         auc_2015_17_media = mean(S_pc))
})) |> mutate(across(where(is.numeric), ~ round(.x, 3)))
# IC95% do AUC interno e da média LOCO por variante (já em 19_auc_ic.csv)
write_csv(sens_tab, file.path(path_tabs, "19_sensibilidade_preditores.csv"))

# ---------------------------------------------------------------------------
# 6) Calibração por subgrupo, limiares de decisão e curva de decisão (modelo completo, teste interno)
# ---------------------------------------------------------------------------
ti <- dc[idx_te, c("pais", "regiao_ssa", "riqueza", "escolaridade", "residencia", "idade_grp")] |>
  mutate(y = y[idx_te], p = res$completo$interno$p, r = 1 - p)
dims <- c("riqueza", "escolaridade", "residencia", "idade_grp")
cal_sub <- bind_rows(lapply(dims, function(dm) {
  ti |> filter(!is.na(.data[[dm]])) |> group_by(dimensao = dm, grupo = as.character(.data[[dm]])) |>
    summarise(n = n(), cobertura_obs = mean(y), cobertura_prev = mean(p), razao_O_E = mean(y) / mean(p),
              cal_int = cal_int(y, p), cal_slope = cal_slope(y, p), brier = brier(y, p), auc = auc_fast(y, p), .groups = "drop")
})) |> mutate(across(where(is.numeric), ~ round(.x, 3)))
write_csv(cal_sub, file.path(path_tabs, "19_calibracao_subgrupo.csv"))

# Limiares: (a) risco >= 0,5; (b) sinalizar os 30% de maior risco (capacidade); (c) sinalizar 50%
medidas_limiar <- function(d, flag) {
  ev <- d$y == 0                                    # evento = não testada
  tibble(n = nrow(d), prop_sinalizada = mean(flag), sens = sum(flag & ev) / sum(ev), espec = sum(!flag & !ev) / sum(!ev),
         vpp = ifelse(sum(flag) > 0, sum(flag & ev) / sum(flag), NA), vpn = ifelse(sum(!flag) > 0, sum(!flag & !ev) / sum(!flag), NA))
}
lim_def <- list(`risco >= 0,5` = quote(ti$r >= 0.5),
                `30% de maior risco` = quote(ti$r >= quantile(ti$r, 0.70)),
                `50% de maior risco` = quote(ti$r >= quantile(ti$r, 0.50)))
lim_tab <- bind_rows(lapply(names(lim_def), function(nm) {
  ti$flag <- eval(lim_def[[nm]])
  bind_rows(
    medidas_limiar(ti, ti$flag) |> mutate(limiar = nm, dimensao = "Geral", grupo = "Todas"),
    bind_rows(lapply(dims, function(dm) {
      ti |> filter(!is.na(.data[[dm]])) |> group_by(dimensao = dm, grupo = as.character(.data[[dm]])) |>
        group_modify(~ medidas_limiar(.x, .x$flag)) |> ungroup() |> mutate(limiar = nm)
    })))
})) |> select(limiar, dimensao, grupo, everything()) |> mutate(across(where(is.numeric), ~ round(.x, 3)))
write_csv(lim_tab, file.path(path_tabs, "19_limiares_subgrupo.csv"))

# Decision curve analysis: evento = não testada; intervenção = ação dirigida às sinalizadas (r >= pt)
nb_modelo <- function(d, pt) { flag <- d$r >= pt; ev <- d$y == 0
  sum(flag & ev) / nrow(d) - sum(flag & !ev) / nrow(d) * pt / (1 - pt) }
nb_todas <- function(d, pt) { ev <- d$y == 0; mean(ev) - mean(!ev) * pt / (1 - pt) }
pts <- seq(0.05, 0.90, by = 0.05)
dca <- bind_rows(
  tibble(grupo = "Todas", pt = pts, modelo = sapply(pts, function(t) nb_modelo(ti, t)), tratar_todas = sapply(pts, function(t) nb_todas(ti, t))),
  bind_rows(lapply(c("Q1 (mais pobre)", "Q5 (mais rico)"), function(g) {
    d <- ti[ti$riqueza %in% g, ]
    tibble(grupo = g, pt = pts, modelo = sapply(pts, function(t) nb_modelo(d, t)), tratar_todas = sapply(pts, function(t) nb_todas(d, t)))
  }))) |> mutate(tratar_nenhuma = 0)
write_csv(dca |> mutate(across(where(is.numeric), ~ round(.x, 4))), file.path(path_tabs, "19_dca.csv"))

dl <- dca |> pivot_longer(c(modelo, tratar_todas, tratar_nenhuma), names_to = "estrategia", values_to = "nb") |>
  mutate(estrategia = recode(estrategia, modelo = "Modelo (gradient boosting)", tratar_todas = "Sinalizar todas", tratar_nenhuma = "Não sinalizar nenhuma"),
         grupo = factor(grupo, levels = c("Todas", "Q1 (mais pobre)", "Q5 (mais rico)")))
g <- ggplot(dl, aes(pt, nb, colour = estrategia, linetype = estrategia)) + geom_line(linewidth = 0.8) +
  facet_wrap(~ grupo, nrow = 1) + coord_cartesian(ylim = c(-0.05, 0.45)) +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "grey40")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  labs(x = "Limiar de risco de não testagem (probabilidade)", y = "Benefício líquido") + tema_pub(11)
ggsave(file.path(path_figs, "19_fig_dca.png"), g, width = 10, height = 3.8, dpi = 300)

# ---------------------------------------------------------------------------
# 7) AUC por região (gradient boosting e random forest) no teste interno
# ---------------------------------------------------------------------------
prf <- ler_parquet(file.path(path_process, "pred_teste_interno.parquet")) |> select(id_mulher, p_ranger, p_xgboost)
tir <- dc[idx_te, c("id_mulher", "regiao_ssa")] |> inner_join(prf, by = "id_mulher") |> mutate(y = y[match(id_mulher, dc$id_mulher)])
reg <- tir |> group_by(regiao = regiao_ssa) |>
  summarise(n = n(), auc_xgboost = auc_fast(y, p_xgboost), auc_ranger = auc_fast(y, p_ranger), .groups = "drop") |>
  mutate(across(starts_with("auc"), ~ round(.x, 3)))
write_csv(reg, file.path(path_tabs, "19_auc_regiao_xgb.csv"))

# ---------------------------------------------------------------------------
# 8) Tamanho de amostra (Riley et al., 2019) para o modelo de desenvolvimento
# ---------------------------------------------------------------------------
Xp <- model.matrix(~ ., data = dc[pred_all])[, -1]; Xp <- Xp[, colSums(Xp[idx_tr, , drop = FALSE]) > 0, drop = FALSE]
n_tr <- length(idx_tr); ev_tr <- sum(y[idx_tr] == 0); phi <- mean(y[idx_tr])
lr_m <- suppressWarnings(glm(y[idx_tr] ~ Xp[idx_tr, ], family = binomial))
lnL0 <- n_tr * (phi * log(phi) + (1 - phi) * log(1 - phi))
R2cs <- 1 - exp(-(as.numeric(lr_m$null.deviance - lr_m$deviance)) / n_tr)
R2max <- 1 - exp(2 * lnL0 / n_tr)
P <- ncol(Xp)
n_crit <- function(S) P / ((S - 1) * log(1 - R2cs / S))
S3 <- R2cs / (R2cs + 0.05 * R2max)
ss_tab <- tibble(
  n_treino = n_tr, eventos_nao_testadas = ev_tr, nao_eventos_testadas = n_tr - ev_tr, parametros_candidatos = P,
  eventos_por_parametro = round(min(ev_tr, n_tr - ev_tr) / P, 1), R2_CoxSnell_aparente = round(R2cs, 3),
  n_min_criterio1_encolhimento_0.9 = ceiling(n_crit(0.9)),
  n_min_criterio2_risco_global_margem_0.05 = ceiling((1.96 / 0.05)^2 * phi * (1 - phi)),
  n_min_criterio3_otimismo_0.05 = ceiling(n_crit(S3)))
ss_tab <- ss_tab |> mutate(n_min_necessario = pmax(n_min_criterio1_encolhimento_0.9, n_min_criterio2_risco_global_margem_0.05, n_min_criterio3_otimismo_0.05))
write_csv(ss_tab, file.path(path_tabs, "19_tamanho_amostra.csv"))

# ---------------------------------------------------------------------------
# 9) Cobertura entre mulheres sem pré-natal (número citado na Discussão)
# ---------------------------------------------------------------------------
sem <- dc |> filter(amostra == "principal", anc_qualquer == "Não") |>
  summarise(n = n(), n_testadas = sum(y), cobertura_nao_ponderada = 100 * mean(y),
            cobertura_w_igual = 100 * weighted.mean(y, w_igual), cobertura_w_dhs = 100 * weighted.mean(y, w_dhs))
com <- dc |> filter(amostra == "principal", anc_qualquer == "Sim") |>
  summarise(n = n(), cobertura_nao_ponderada = 100 * mean(y), cobertura_w_igual = 100 * weighted.mean(y, w_igual), cobertura_w_dhs = 100 * weighted.mean(y, w_dhs))
write_csv(bind_rows(mutate(sem, grupo = "Sem pré-natal"), mutate(com, grupo = "Com pré-natal")) |> mutate(across(where(is.numeric), ~ round(.x, 2))),
          file.path(path_tabs, "19_testagem_sem_cpn.csv"))

cat("\nEtapa 19 concluída.\n")
cat("-- AUC agregado vs média por país:\n"); print(as.data.frame(agr))
cat("-- Brier escalonado:\n"); print(as.data.frame(bs_tab))
cat("-- Recalibração (CV por PSU):\n"); print(as.data.frame(rc))
cat("-- Sensibilidade a preditores:\n"); print(as.data.frame(sens_tab))
cat("-- Calibração por subgrupo:\n"); print(as.data.frame(cal_sub))
cat("-- Limiares (geral e riqueza):\n"); print(as.data.frame(lim_tab |> filter(dimensao %in% c("Geral", "riqueza"))))
cat("-- Tamanho de amostra:\n"); print(as.data.frame(t(ss_tab)))
cat("-- AUC por região:\n"); print(as.data.frame(reg))
cat("-- Sem CPN:\n"); print(as.data.frame(sem))
