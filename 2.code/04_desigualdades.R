# ============================================================================
# P06 - Testagem de HIV no pré-natal na África Subsaariana
# Etapa 04 - Desigualdades na cobertura de testagem (medidas de iniquidade)
# Autor: Manuel Mahoche
# Entrada: 1.data/process/ds_clean.parquet
# Saída:   3.output/tabs/04_sii_rii_concentracao.csv  SII, RII e índice de concentração
#          3.output/tabs/04_gaps_extremos.csv          gaps absolutos e razões (categorias extremas)
#
# Desenho: igual à Etapa 03 (svydesign com psu_id/estrato_id).
#   - POR PAÍS: peso w_norm (representativo dentro do país; rank e medidas por país).
#   - POOLED (SSA): peso w_igual (cada país contribui igualmente); rank sobre a amostra agrupada.
#
# Medidas (padrão WHO Health Equity Monitor):
#   - SII (Slope Index of Inequality): dif. absoluta de cobertura entre topo e base do
#     gradiente socioeconómico, via regressão do desfecho no rank fracionário (ridit). Em pp.
#   - RII (Relative Index of Inequality): razão topo/base, via regressão de Poisson no rank.
#   - Índice de concentração (riqueza): regressão conveniente de Kakwani (SE do desenho).
#     Positivo = cobertura concentrada nos mais ricos (pró-rico).
#   - Gaps extremos: dif. (pp) e razão entre as categorias extremas, via svyglm com a
#     dimensão como fator (ref = categoria de base).
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(survey); library(here)
})
source(here("2.code", "R", "config.R"))
source(here("2.code", "R", "funcoes.R"))
options(survey.lonely.psu = "adjust")

dc <- ler_parquet(file.path(path_process, "ds_clean.parquet")) |>
  filter(amostra == "principal") |>
  mutate(desfecho_num = as.numeric(hiv_anc_resultado))

# Rank fracionário (ridit) ponderado de uma variável ordinal
ridit <- function(gi, w) {
  ws <- tapply(w, gi, sum); ws <- ws[order(as.numeric(names(ws)))]
  tot <- sum(ws); cumprev <- cumsum(c(0, ws))[seq_along(ws)]; names(cumprev) <- names(ws)
  r <- (cumprev + ws / 2) / tot
  r[as.character(gi)]
}
# Extrai estimativa + IC95% (Wald) de um coeficiente, com transformação opcional.
# IC calculado a partir de summary()$coefficients (Estimate ± t * SE, graus de liberdade do
# desenho). Funciona para gaussiana e quasipoisson sem depender do despacho de confint().
getco <- function(m, nm, trans = function(x) x) {
  if (is.null(m)) return(c(est = NA, lo = NA, hi = NA))
  co <- tryCatch(summary(m)$coefficients, error = function(e) NULL)
  if (is.null(co) || !(nm %in% rownames(co))) return(c(est = NA, lo = NA, hi = NA))
  b <- co[nm, 1]; se <- co[nm, 2]
  df <- tryCatch(stats::df.residual(m), error = function(e) Inf)
  if (!is.finite(df) || df <= 0) df <- Inf
  q <- stats::qt(0.975, df)
  stats::setNames(trans(c(b, b - q * se, b + q * se)), c("est", "lo", "hi"))
}
f1 <- function(x) sprintf("%.1f", x)
f2 <- function(x) sprintf("%.2f", x)
fic <- function(v, f) ifelse(is.na(v["est"]), "—", sprintf("%s (%s–%s)", f(v["est"]), f(v["lo"]), f(v["hi"])))

# Dimensões socioeconómicas ordenadas (SII/RII/concentração) e extremos para os gaps
dims_grad <- list(
  riqueza      = list(top = "riquezaQ5 (mais rico)",            contraste = "Q5 vs Q1",      conc = TRUE),
  escolaridade = list(top = "escolaridadeSecundária ou mais",   contraste = "Secundária+ vs Nenhuma", conc = FALSE)
)
# Gaps também para residência (não ordinal; só dif/razão entre categorias)
dims_gap_extra <- list(
  residencia = list(top = "residenciaRural", contraste = "Rural vs Urbana")
)

# ---------------------------------------------------------------------------
# Núcleo: medidas para um âmbito (país ou pooled), dado o peso a usar
# ---------------------------------------------------------------------------
medidas_ambito <- function(data, wcol, ambito) {
  sr_rows <- list(); gap_rows <- list()

  # --- SII / RII / concentração (dimensões graduadas) ---
  for (dim in names(dims_grad)) {
    cfg <- dims_grad[[dim]]
    d <- data |> filter(!is.na(.data[[dim]]), !is.na(desfecho_num))
    if (nrow(d) < 50 || dplyr::n_distinct(d[[dim]]) < 2) next
    d$gi    <- as.integer(d[[dim]])
    d$rankv <- ridit(d$gi, d[[wcol]])
    w <- d[[wcol]]; y <- d$desfecho_num
    mu <- weighted.mean(y, w)
    s2r <- sum(w * (d$rankv - weighted.mean(d$rankv, w))^2) / sum(w)
    d$lhs_ci <- 2 * s2r * (y / mu)
    des <- svydesign(ids = ~psu_id, strata = ~estrato_id,
                     weights = as.formula(paste0("~", wcol)), data = d, nest = TRUE)
    m_sii <- tryCatch(svyglm(desfecho_num ~ rankv, design = des), error = function(e) NULL)
    m_rii <- tryCatch(svyglm(desfecho_num ~ rankv, design = des, family = quasipoisson()), error = function(e) NULL)
    # Nomes de variável distintos dos nomes de coluna: dentro de tibble() uma coluna já
    # criada ofusca a variável homónima (era a causa de rii_ic sair "(NA–NA)").
    sii_v <- getco(m_sii, "rankv", function(x) 100 * x)
    rii_v <- getco(m_rii, "rankv", function(x) exp(x))
    ci_v  <- c(est = NA, lo = NA, hi = NA)
    if (isTRUE(cfg$conc)) {
      m_ci <- tryCatch(svyglm(lhs_ci ~ rankv, design = des), error = function(e) NULL)
      ci_v <- getco(m_ci, "rankv")
    }
    sr_rows[[dim]] <- tibble(
      ambito = ambito, dimensao = dim, n = nrow(d),
      sii_pp = round(sii_v["est"], 1), sii_ic = fic(sii_v, f1),
      rii = round(rii_v["est"], 2), rii_ic = fic(rii_v, f2),
      cindex = round(ci_v["est"], 3),
      cindex_ic = if (isTRUE(cfg$conc)) fic(ci_v, function(x) sprintf("%.3f", x)) else "")
  }

  # --- Gaps extremos (dif pp e razão) para todas as dimensões ---
  for (dim in names(c(dims_grad, dims_gap_extra))) {
    cfg <- c(dims_grad, dims_gap_extra)[[dim]]
    d <- data |> filter(!is.na(.data[[dim]]), !is.na(desfecho_num))
    if (nrow(d) < 50 || dplyr::n_distinct(d[[dim]]) < 2) next
    des <- svydesign(ids = ~psu_id, strata = ~estrato_id,
                     weights = as.formula(paste0("~", wcol)), data = d, nest = TRUE)
    m_dif <- tryCatch(svyglm(as.formula(paste0("desfecho_num ~ `", dim, "`")), design = des),
                      error = function(e) NULL)
    m_raz <- tryCatch(svyglm(as.formula(paste0("desfecho_num ~ `", dim, "`")), design = des,
                             family = quasipoisson()), error = function(e) NULL)
    dif <- getco(m_dif, cfg$top, function(x) 100 * x)
    raz <- getco(m_raz, cfg$top, function(x) exp(x))
    gap_rows[[dim]] <- tibble(
      ambito = ambito, dimensao = dim, contraste = cfg$contraste, n = nrow(d),
      dif_pp = round(dif["est"], 1), dif_ic = fic(dif, f1),
      razao = round(raz["est"], 2), razao_ic = fic(raz, f2))
  }
  list(sr = bind_rows(sr_rows), gap = bind_rows(gap_rows))
}

# ---------------------------------------------------------------------------
# Corre por país (w_norm) e pooled (w_igual)
# ---------------------------------------------------------------------------
paises <- dc |> distinct(pais, nome_pais) |> arrange(pais)
res_pais <- lapply(paises$pais, function(p) {
  r <- medidas_ambito(filter(dc, pais == p), "w_norm",
                      paises$nome_pais[paises$pais == p])
  r
})
res_pool <- medidas_ambito(dc, "w_igual", "SSA (pooled, w_igual)")

sii_rii <- bind_rows(res_pool$sr, bind_rows(lapply(res_pais, `[[`, "sr"))) |>
  mutate(dimensao = recode(dimensao, riqueza = "Riqueza", escolaridade = "Escolaridade"),
         # SII e RII são extrapolações do gradiente para os extremos do rank (0 e 1);
         # em países de cobertura baixa e gradiente íngreme sobre-extrapolam (SII > 100 pp,
         # RII explode porque a base extrapolada tende a zero). Marcado para leitura cautelosa.
         # O índice de concentração e o gap Q5–Q1 (04_gaps_extremos.csv) são as medidas robustas.
         nota = dplyr::if_else(abs(sii_pp) > 95 | rii > 5,
                               "extrapolação instável (cobertura baixa) — usar índice de concentração/gap", ""))
gaps <- bind_rows(res_pool$gap, bind_rows(lapply(res_pais, `[[`, "gap"))) |>
  mutate(dimensao = recode(dimensao, riqueza = "Riqueza", escolaridade = "Escolaridade",
                           residencia = "Residência"))

write_csv(sii_rii, file.path(path_tabs, "04_sii_rii_concentracao.csv"))
write_csv(gaps,    file.path(path_tabs, "04_gaps_extremos.csv"))

# Resumo no ecrã: concentração por riqueza, do mais pró-rico ao menos
resumo <- sii_rii |> filter(dimensao == "Riqueza", ambito != "SSA (pooled, w_igual)") |>
  arrange(desc(cindex))
cat("Etapa 04 concluída.\n")
cat(sprintf("  Pooled (w_igual) riqueza: SII=%.1f pp, RII=%.2f, índice de concentração=%.3f\n",
            sii_rii$sii_pp[sii_rii$ambito == "SSA (pooled, w_igual)" & sii_rii$dimensao == "Riqueza"],
            sii_rii$rii[sii_rii$ambito == "SSA (pooled, w_igual)" & sii_rii$dimensao == "Riqueza"],
            sii_rii$cindex[sii_rii$ambito == "SSA (pooled, w_igual)" & sii_rii$dimensao == "Riqueza"]))
cat("  Países mais pró-rico (índice de concentração por riqueza):\n")
print(utils::head(resumo[, c("ambito", "sii_pp", "rii", "cindex")], 5))
cat("  Tabelas: 04_sii_rii_concentracao.csv, 04_gaps_extremos.csv\n")
