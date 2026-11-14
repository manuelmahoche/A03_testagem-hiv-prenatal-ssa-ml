# ============================================================================
# P06 - Testagem de HIV no pré-natal na África Subsaariana
# Etapa 18 - Desigualdade: índice de concentração corrigido (Wagstaff, Erreygers)
# Autor: Manuel Mahoche
# Entrada: 1.data/process/ds_clean.parquet
# Saída:   3.output/tabs/18_desigualdade_corrigida.csv   CI bruto/Wagstaff/Erreygers, SII, RII por país
#          3.output/tabs/18_desigualdade_vs_cobertura.csv correlação de Spearman medida x cobertura
#          3.output/figs/18_fig_desigualdade_cobertura.png
#
# Motivo (revisão interna): para desfechos binários o índice de concentração bruto C depende
# da cobertura média mu (limite superior = 1 - mu), o que o infla nos países de baixa
# cobertura. Correções para desfechos limitados (O'Donnell et al. 2008, cap. 8):
#   Wagstaff:  W = C / (1 - mu)
#   Erreygers: E = 4 * mu * C         (preferida; satisfaz o critério de espelho)
# SII/RII mantêm-se como em 04 para comparação.
# Posição socioeconómica: ridit do quintil de riqueza do DHS (v190), calculado DENTRO de
# cada país (o índice de riqueza do DHS é específico de cada país). No agregado (peso
# w_igual) usa-se o ridit calculado dentro do país, e não o ridit sobre os quintis empilhados.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(survey); library(here); library(ggplot2)
})
source(here("2.code", "R", "config.R"))
source(here("2.code", "R", "funcoes.R"))
source(here("2.code", "R", "tema_figuras.R"))
options(survey.lonely.psu = "adjust")

dc <- ler_parquet(file.path(path_process, "ds_clean.parquet")) |>
  filter(amostra == "principal", !is.na(riqueza), !is.na(hiv_anc_resultado)) |>
  mutate(y = as.numeric(hiv_anc_resultado), gi = as.integer(riqueza))

ridit <- function(gi, w) {
  ws <- tapply(w, gi, sum); ws <- ws[order(as.numeric(names(ws)))]
  tot <- sum(ws); cumprev <- cumsum(c(0, ws))[seq_along(ws)]; names(cumprev) <- names(ws)
  r <- (cumprev + ws / 2) / tot
  unname(r[as.character(gi)])
}
# ridit dentro de cada país (peso w_norm = representativo no país)
dc <- dc |> group_by(pais) |> mutate(rank_pais = ridit(gi, w_norm)) |> ungroup()

# Estimativa + IC95% (Wald, gl do desenho) de um coeficiente
coef_ic <- function(m, nm) {
  co <- summary(m)$coefficients; b <- co[nm, 1]; se <- co[nm, 2]
  df <- stats::df.residual(m); if (!is.finite(df) || df <= 0) df <- Inf
  q <- qt(0.975, df); c(est = b, se = se, lo = b - q * se, hi = b + q * se)
}

medidas <- function(d, wcol, rankcol, rotulo) {
  w <- d[[wcol]]; y <- d$y; r <- d[[rankcol]]
  mu <- weighted.mean(y, w)
  s2r <- sum(w * (r - weighted.mean(r, w))^2) / sum(w)
  d$rk <- r; d$lhs <- 2 * s2r * (y / mu)
  des <- svydesign(ids = ~psu_id, strata = ~estrato_id, weights = as.formula(paste0("~", wcol)),
                   data = d, nest = TRUE)
  ci  <- coef_ic(svyglm(lhs ~ rk, design = des), "rk")
  sii <- coef_ic(svyglm(y ~ rk, design = des), "rk") * 100
  rii_v <- tryCatch({
    m <- svyglm(y ~ rk, design = des, family = quasipoisson()); exp(coef_ic(m, "rk")[c("est", "lo", "hi")])
  }, error = function(e) c(est = NA, lo = NA, hi = NA))
  tibble(ambito = rotulo, n = nrow(d), cobertura_pct = round(100 * mu, 1),
         ci_bruto = ci["est"], ci_lo = ci["lo"], ci_hi = ci["hi"],
         wag = ci["est"] / (1 - mu), wag_lo = ci["lo"] / (1 - mu), wag_hi = ci["hi"] / (1 - mu),
         err = 4 * mu * ci["est"], err_lo = 4 * mu * ci["lo"], err_hi = 4 * mu * ci["hi"],
         sii_pp = sii["est"], sii_lo = sii["lo"], sii_hi = sii["hi"],
         rii = rii_v["est"], rii_lo = rii_v["lo"], rii_hi = rii_v["hi"])
}

paises <- dc |> distinct(pais, nome_pais, estrato) |> arrange(pais)
res_pais <- bind_rows(lapply(seq_len(nrow(paises)), function(i) {
  d <- filter(dc, pais == paises$pais[i])
  medidas(d, "w_norm", "rank_pais", paises$nome_pais[i]) |> mutate(pais = paises$pais[i], estrato = paises$estrato[i])
}))
# Agregado: peso w_igual, ridit dentro do país
res_pool <- medidas(dc, "w_igual", "rank_pais", "SSA (agregado, w_igual; ridit dentro do país)") |>
  mutate(pais = "POOL", estrato = NA_character_)
# Agregado alternativo (ridit sobre os quintis empilhados, como na etapa 04)
dc$rank_emp <- ridit(dc$gi, dc$w_igual)
res_pool_emp <- medidas(dc, "w_igual", "rank_emp", "SSA (agregado, w_igual; ridit sobre quintis empilhados)") |>
  mutate(pais = "POOL_EMP", estrato = NA_character_)

out <- bind_rows(res_pool, res_pool_emp, res_pais) |>
  mutate(across(c(ci_bruto, ci_lo, ci_hi, wag, wag_lo, wag_hi, err, err_lo, err_hi), ~ round(.x, 3)),
         across(c(sii_pp, sii_lo, sii_hi), ~ round(.x, 1)),
         across(c(rii, rii_lo, rii_hi), ~ round(.x, 2))) |>
  select(pais, ambito, estrato, n, cobertura_pct, everything())
write_csv(out, file.path(path_tabs, "18_desigualdade_corrigida.csv"))

# Padrão desigualdade x cobertura (23 países): Spearman
p <- filter(out, !pais %in% c("POOL", "POOL_EMP"))
sp <- bind_rows(lapply(c(ci_bruto = "ci_bruto", wagstaff = "wag", erreygers = "err", sii_pp = "sii_pp"), function(v) {
  t <- suppressWarnings(cor.test(p$cobertura_pct, p[[v]], method = "spearman", exact = FALSE))
  tibble(rho = unname(t$estimate), p = t$p.value)
}), .id = "medida") |> mutate(rho = round(rho, 2), p = signif(p, 2))
# Erreygers por tercil de cobertura
p <- p |> mutate(tercil = ntile(cobertura_pct, 3))
terc <- p |> group_by(tercil) |> summarise(n_paises = n(), cob_min = min(cobertura_pct), cob_max = max(cobertura_pct),
                                           ci_bruto_med = round(median(ci_bruto), 3), err_med = round(median(err), 3),
                                           sii_med = round(median(sii_pp), 1), .groups = "drop")
write_csv(sp, file.path(path_tabs, "18_desigualdade_vs_cobertura.csv"))
write_csv(terc, file.path(path_tabs, "18_desigualdade_tercis_cobertura.csv"))

# Figura: três painéis (bruto, Erreygers, SII) contra a cobertura
lg <- p |> select(nome = ambito, cobertura_pct, `Índice de concentração bruto` = ci_bruto,
                  `Índice de Erreygers` = err, `SII (pp)` = sii_pp) |>
  pivot_longer(-c(nome, cobertura_pct), names_to = "medida", values_to = "valor") |>
  mutate(medida = factor(medida, levels = c("Índice de concentração bruto", "Índice de Erreygers", "SII (pp)")))
g <- ggplot(lg, aes(cobertura_pct, valor)) +
  geom_point(size = 2, colour = "#0072B2") +
  geom_smooth(method = "loess", se = FALSE, colour = "#D55E00", linewidth = 0.7, formula = y ~ x) +
  facet_wrap(~ medida, scales = "free_y", nrow = 1) +
  labs(x = "Cobertura de testagem de HIV no pré-natal (%)", y = NULL) +
  tema_pub(11)
ggsave(file.path(path_figs, "18_fig_desigualdade_cobertura.png"), g, width = 10, height = 3.6, dpi = 300)

cat("Etapa 18 concluída.\n")
cat("  Agregado (ridit dentro do país vs empilhado):\n")
print(as.data.frame(out[out$pais %in% c("POOL", "POOL_EMP"), c("ambito", "cobertura_pct", "ci_bruto", "wag", "err", "sii_pp", "rii")]))
cat("  Spearman medida x cobertura (23 países):\n"); print(as.data.frame(sp))
cat("  Por tercil de cobertura:\n"); print(as.data.frame(terc))
print(as.data.frame(p[order(p$cobertura_pct), c("ambito", "cobertura_pct", "ci_bruto", "wag", "err", "sii_pp")]))
