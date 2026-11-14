# ============================================================================
# P06 - Testagem de HIV no pré-natal na África Subsaariana
# Etapa 11 - Sensibilidade a dados em falta e ao desenho (análogo ao ESM7 da referência)
# Autor: Manuel Mahoche
# Entrada: ds_clean.parquet, particoes.parquet, 06_modelos_finais.rds
# Saída:   3.output/tabs/11_caso_completo_vs_ausente.csv  características: completos vs ausentes
#          3.output/tabs/11_padrao_ausencia.csv           ausência por variável e co-ausência
#          3.output/tabs/11_sensibilidade_na.csv          desempenho por estratégia de NA
#          3.output/tabs/11_ponderado_vs_nao.csv          treino ponderado vs não (AUC + ranks)
#
# Replica, adaptado à amostra multi-país (Estrato A), as análises do ESM7:
#   - Tabela Sa: características (não ponderadas) de casos completos vs com ≥1 preditor ausente.
#   - Tabela Sb: padrão de ausência (quais variáveis e co-ocorrência).
#   - Tabela Sc (estratégia de NA): desempenho no teste do xgboost sob 3 estratégias —
#       primária "(sem dado)" como categoria; caso-completo; imputação pela moda.
#     (Imputação múltipla do original aproximada por imputação única pela moda, por viabilidade.)
#   - Tabelas Sd/Se: treino não ponderado vs ponderado (w_norm) — AUC e estabilidade do ranking.
#       SMOTE do original não é replicado: o desfecho é ~42% (quase equilibrado), onde o
#       reequilíbrio de classes acrescenta pouco; documentado aqui.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(here); library(xgboost); library(pROC)
})
source(here("2.code", "R", "config.R"))
source(here("2.code", "R", "funcoes.R"))
set.seed(semente)

pred <- setdiff(vars_preditores, "regiao_id")
dc <- ler_parquet(file.path(path_process, "ds_clean.parquet")) |>
  filter(amostra == "principal", estrato == "A") |>
  select(id_mulher, pais, all_of(pred), w_norm, y = hiv_anc_resultado)
pa <- ler_parquet(file.path(path_process, "particoes.parquet")) |> select(id_mulher, particao)
dc <- inner_join(dc, pa, by = "id_mulher")
hp <- readRDS(file.path("3.output", "mod", "06_modelos_finais.rds"))$hp_xgb

dc$n_miss <- rowSums(is.na(dc[pred]))
dc$grupo_na <- ifelse(dc$n_miss > 0, "Com ausência", "Caso completo")

# --- Tabela Sa: características (não ponderadas) completos vs com ausência ---
sa <- lapply(pred, function(v) {
  d <- dc |> filter(!is.na(.data[[v]]))
  tb <- d |> count(grupo_na, cat = .data[[v]]) |>
    group_by(grupo_na) |> mutate(pct = round(100 * n / sum(n), 1)) |> ungroup()
  pval <- tryCatch(suppressWarnings(chisq.test(table(d$grupo_na, d[[v]]))$p.value), error = function(e) NA_real_)
  tb |> mutate(variavel = v, p_valor = ifelse(cat == first(cat), format.pval(pval, digits = 2, eps = 1e-4), "")) |>
    select(variavel, categoria = cat, grupo_na, n, pct, p_valor)
}) |> bind_rows()
write_csv(sa, file.path(path_tabs, "11_caso_completo_vs_ausente.csv"))

# --- Tabela Sb: padrão de ausência ---
sb_var <- tibble(variavel = pred,
                 n_ausente = sapply(pred, function(v) sum(is.na(dc[[v]]))),
                 pct_ausente = round(100 * sapply(pred, function(v) mean(is.na(dc[[v]]))), 2)) |>
  arrange(desc(n_ausente))
sb_pad <- dc |> count(n_vars_ausentes = n_miss, name = "n_mulheres") |>
  mutate(pct = round(100 * n_mulheres / sum(n_mulheres), 2)) |> arrange(n_vars_ausentes)
write_csv(bind_rows(
  sb_var |> mutate(tipo = "por_variavel") |> rename(chave = variavel, n = n_ausente, pct = pct_ausente),
  sb_pad |> mutate(tipo = "n_vars_ausentes_por_mulher", chave = as.character(n_vars_ausentes)) |>
    rename(n = n_mulheres) |> select(tipo, chave, n, pct)),
  file.path(path_tabs, "11_padrao_ausencia.csv"))

# --- Funções de modelagem (xgboost com hp da Etapa 06) ---
fill_sem_dado <- function(d) { for (v in pred) { x <- as.factor(d[[v]]); if (anyNA(x)) { x <- addNA(x); levels(x)[is.na(levels(x))] <- "(sem dado)" }; d[[v]] <- x }; d }
fill_moda <- function(d) { for (v in pred) { x <- as.factor(d[[v]]); if (anyNA(x)) { m <- names(which.max(table(x))); x[is.na(x)] <- m; x <- droplevels(x) }; d[[v]] <- x }; d }

treina_avalia <- function(train, test, w = NULL) {
  lv <- lapply(train[pred], levels)
  Xtr <- model.matrix(~ . - 1, data = train[pred])
  Xte <- model.matrix(~ . - 1, data = test[pred])[, colnames(Xtr), drop = FALSE]
  dtr <- if (is.null(w)) xgb.DMatrix(Xtr, label = train$y) else xgb.DMatrix(Xtr, label = train$y, weight = w)
  m <- xgb.train(list(objective = "binary:logistic", eval_metric = "auc", max_depth = hp$max_depth,
                      eta = hp$eta, subsample = 0.8, colsample_bytree = 0.8, nthread = 0),
                 dtr, nrounds = hp$nrounds, verbose = 0)
  p <- as.numeric(predict(m, Xte))
  list(auc = as.numeric(pROC::auc(pROC::roc(test$y, p, quiet = TRUE, direction = "<"))),
       brier = mean((p - test$y)^2), model = m, cols = colnames(Xtr))
}

# --- Tabela Sc: estratégia de dados em falta ---
estrategias <- list(
  `Primária: "(sem dado)"` = list(prep = fill_sem_dado, completo = FALSE),
  `Caso completo`          = list(prep = identity,      completo = TRUE),
  `Imputação pela moda`    = list(prep = fill_moda,     completo = FALSE))
sc <- lapply(names(estrategias), function(nm) {
  cfg <- estrategias[[nm]]
  d <- dc; if (cfg$completo) d <- d |> filter(n_miss == 0)
  d <- cfg$prep(d)
  tr <- d |> filter(particao == "treino"); te <- d |> filter(particao == "teste")
  r <- treina_avalia(tr, te)
  tibble(estrategia = nm, n_treino = nrow(tr), n_teste = nrow(te),
         auc = round(r$auc, 3), brier = round(r$brier, 3))
}) |> bind_rows()
write_csv(sc, file.path(path_tabs, "11_sensibilidade_na.csv"))

# --- Tabelas Sd/Se: ponderado (w_norm) vs não ponderado, estratégia primária ---
dps <- fill_sem_dado(dc); trp <- dps |> filter(particao == "treino"); tep <- dps |> filter(particao == "teste")
r_unw <- treina_avalia(trp, tep, w = NULL)
r_wei <- treina_avalia(trp, tep, w = trp$w_norm)
imp_rank <- function(m, cols) {
  imp <- xgb.importance(model = m); g <- setNames(rep(0, length(cols)), cols)
  g[imp$Feature] <- imp$Gain
  # agrega por variável (prefixo mais longo)
  v <- vapply(names(g), function(cn) { h <- pred[startsWith(cn, pred)]; h[which.max(nchar(h))] }, character(1))
  tapply(g, v, sum)
}
gu <- imp_rank(r_unw$model, r_unw$cols); gw <- imp_rank(r_wei$model, r_wei$cols)
ru <- rank(-gu); rw <- rank(-gw[names(gu)])
se <- tibble(variavel = names(gu),
             rank_nao_ponderado = as.integer(ru), rank_ponderado = as.integer(rw[names(gu)]),
             dif_rank = rank_ponderado - rank_nao_ponderado) |> arrange(rank_nao_ponderado)
write_csv(
  bind_rows(
    tibble(metrica = "AUC teste", nao_ponderado = round(r_unw$auc, 3), ponderado = round(r_wei$auc, 3)),
    tibble(metrica = "Brier teste", nao_ponderado = round(r_unw$brier, 3), ponderado = round(r_wei$brier, 3))),
  file.path(path_tabs, "11_ponderado_vs_nao.csv"))
write_csv(se, file.path(path_tabs, "11_ranks_ponderado_vs_nao.csv"))

cat("Etapa 11 concluída.\n")
cat("  Estratégia de dados em falta (AUC teste):\n"); print(as.data.frame(sc))
cat(sprintf("  Treino: AUC não ponderado %.3f vs ponderado %.3f.\n", r_unw$auc, r_wei$auc))
cat(sprintf("  Mulheres com >=1 preditor ausente: %d (%.1f%%).\n",
            sum(dc$n_miss > 0), 100 * mean(dc$n_miss > 0)))
cat("  Tabelas: 11_caso_completo_vs_ausente, 11_padrao_ausencia, 11_sensibilidade_na, 11_ponderado_vs_nao, 11_ranks_ponderado_vs_nao\n")
