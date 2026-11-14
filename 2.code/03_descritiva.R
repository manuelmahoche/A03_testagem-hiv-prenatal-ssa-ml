# ============================================================================
# P06 - Testagem de HIV no pré-natal na África Subsaariana
# Etapa 03 - Descritiva com desenho amostral (survey)
# Autor: Manuel Mahoche
# Entrada: 1.data/process/ds_clean.parquet
# Saída:   3.output/tabs/03_resumo_amostra.csv        panorama da amostra + prevalência pooled
#          3.output/tabs/03_desfecho_por_pais.csv      prevalência ponderada por país (IC95%)
#          3.output/tabs/03_tabela1_desfecho.csv       Tabela 1: distribuição por desfecho + p + cobertura
#
# Decisões de desenho (fixadas com o autor):
#   - Escopo: amostra PRINCIPAL (23 países, estratos A+B).
#   - Ponderação pooled: w_igual (cada país contribui igualmente para a média SSA).
#   - Tabela 1: coluna geral + por desfecho (testada vs não) + p baseado no desenho.
#
# Desenho: svydesign(ids = psu_id, strata = estrato_id, weights = w_igual, nest = TRUE).
#   psu_id/estrato_id têm prefixo do país, pelo que são únicos entre inquéritos.
#   w_igual já é constante dentro de cada país, logo as estimativas POR PAÍS coincidem
#   com o peso DHS representativo (bate com o STATcompiler); só o pooled fica equilibrado.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(forcats)
  library(survey); library(srvyr); library(here)
})
source(here("2.code", "R", "config.R"))
source(here("2.code", "R", "funcoes.R"))

options(survey.lonely.psu = "adjust")   # estratos com um só PSU: variância ajustada, não erro

dc <- ler_parquet(file.path(path_process, "ds_clean.parquet")) |>
  filter(amostra == "principal") |>
  mutate(desfecho_num = as.numeric(hiv_anc_resultado),
         desfecho_f   = factor(hiv_anc_resultado, levels = c(0, 1),
                               labels = c("Não testada", "Testada")))

verificar(all(!is.na(dc$desfecho_num)), "desfecho com NA na amostra principal")

# Desenho survey (base e srvyr)
des <- svydesign(ids = ~psu_id, strata = ~estrato_id, weights = ~w_igual,
                 data = dc, nest = TRUE)
sd  <- as_survey_design(des)

# Variáveis da Tabela 1, com bloco (segue o dicionário da Etapa 02).
# regiao_id fica de fora: é identificador de desenho com centenas de níveis (região × país),
# sem leitura num quadro agrupado (entra depois como ajuste/efeito fixo nos modelos).
vars_tab <- setdiff(c(vars_preditores, "decisao_saude", vars_pos_exposicao), "regiao_id")
bloco_de <- c(
  setNames(rep("Sociodemográfico",  length(vars_bloco_sociodem)), vars_bloco_sociodem),
  setNames(rep("Informação",        length(vars_bloco_info)),     vars_bloco_info),
  setNames(rep("Reprodutivo",       length(vars_bloco_reprod)),   vars_bloco_reprod),
  setNames(rep("Acesso",            length(vars_bloco_acesso)),   vars_bloco_acesso),
  "decisao_saude" = "Acesso",
  setNames(rep("Pós-exposição (serviço)", length(vars_pos_exposicao)), vars_pos_exposicao)
)
rotulo_var <- c(
  idade_grp = "Faixa etária", uniao = "Estado conjugal", residencia = "Residência",
  regiao_id = "Região", chefe_sexo = "Sexo do chefe do agregado", riqueza = "Quintil de riqueza",
  escolaridade = "Escolaridade", ocupacao = "Ocupação", educ_parceiro = "Escolaridade do parceiro",
  ocup_parceiro = "Ocupação do parceiro", media_semanal = "Exposição semanal a média",
  internet_12m = "Usou internet (12 meses)", conhece_ptv = "Conhece PTV",
  estigma = "Atitude de estigma", paridade = "Paridade", idade_1parto = "Idade ao 1º parto",
  intencao_gravidez = "Intenção da gravidez", desejo_filhos = "Desejo de mais filhos",
  distancia_problema = "Distância como problema", seguro_saude = "Seguro de saúde",
  decisao_saude = "Decisão sobre a própria saúde", anc_qualquer = "Fez algum pré-natal",
  anc_4mais = "4+ consultas de pré-natal", anc_inicio_trim = "Início do pré-natal",
  parto_instituicao = "Parto em instituição"
)
pct1  <- function(x) round(100 * x, 1)
ic_fmt <- function(est, lo, hi) ifelse(is.na(est), "—",
                                       sprintf("%.1f (%.1f–%.1f)", pct1(est), pct1(lo), pct1(hi)))

# ---------------------------------------------------------------------------
# 1) Resumo da amostra + prevalência pooled (w_igual) por estrato
# ---------------------------------------------------------------------------
prev_global <- sd |> summarise(p = survey_mean(desfecho_num, vartype = "ci"))
prev_estr   <- sd |> group_by(estrato) |>
  summarise(p = survey_mean(desfecho_num, vartype = "ci"), n = unweighted(n()))

resumo <- bind_rows(
  tibble(ambito = "Amostra principal (A+B)",
         n_mulheres = nrow(dc), n_paises = n_distinct(dc$pais),
         n_psu = n_distinct(dc$psu_id), n_estratos = n_distinct(dc$estrato_id),
         prev_pct = pct1(prev_global$p),
         prev_ic = sprintf("%.1f–%.1f", pct1(prev_global$p_low), pct1(prev_global$p_upp))),
  prev_estr |> transmute(
    ambito = paste0("Estrato ", estrato), n_mulheres = n,
    n_paises = NA_integer_, n_psu = NA_integer_, n_estratos = NA_integer_,
    prev_pct = pct1(p), prev_ic = sprintf("%.1f–%.1f", pct1(p_low), pct1(p_upp)))
)
write_csv(resumo, file.path(path_tabs, "03_resumo_amostra.csv"))

# ---------------------------------------------------------------------------
# 2) Prevalência do desfecho por país (IC95%), ordenada
# ---------------------------------------------------------------------------
por_pais <- sd |> group_by(pais, nome_pais, regiao_ssa, estrato) |>
  summarise(prev = survey_mean(desfecho_num, vartype = "ci"),
            n = unweighted(n()), .groups = "drop") |>
  transmute(pais, nome_pais, regiao_ssa, estrato, n,
            prev_pct = pct1(prev),
            prev_ic95 = sprintf("%.1f–%.1f", pct1(prev_low), pct1(prev_upp))) |>
  arrange(desc(prev_pct))
write_csv(por_pais, file.path(path_tabs, "03_desfecho_por_pais.csv"))

# ---------------------------------------------------------------------------
# 3) Tabela 1: distribuição geral e por desfecho (coluna %), p e cobertura (linha %)
# ---------------------------------------------------------------------------
linha_cat <- function(v) {
  sdv <- sd |> filter(!is.na(.data[[v]])) |> mutate(cat = .data[[v]])
  # proportion = FALSE: IC linearizado (rápido); suficiente para quadro descritivo.
  prop_por <- function(design) design |> group_by(cat) |>
    summarise(p = survey_prop(vartype = "ci", proportion = FALSE), .groups = "drop")
  g  <- prop_por(sdv) |> rename(geral_p = p, geral_lo = p_low, geral_hi = p_upp)
  dn <- prop_por(sdv |> filter(desfecho_f == "Não testada")) |>
    rename(nt_p = p, nt_lo = p_low, nt_hi = p_upp)
  dt <- prop_por(sdv |> filter(desfecho_f == "Testada")) |>
    rename(t_p = p, t_lo = p_low, t_hi = p_upp)
  # cobertura: % testada dentro de cada categoria (linha %)
  cob <- sdv |> group_by(cat) |>
    summarise(c = survey_mean(desfecho_num, vartype = "ci"), .groups = "drop") |>
    rename(cob_p = c, cob_lo = c_low, cob_hi = c_upp)
  nn <- dc |> filter(!is.na(.data[[v]])) |> count(cat = .data[[v]], name = "n") |>
    mutate(cat = as.character(cat))
  pval <- tryCatch(
    svychisq(as.formula(paste0("~`", v, "` + desfecho_f")), des, statistic = "F")$p.value,
    error = function(e) NA_real_)
  g |> mutate(cat = as.character(cat)) |>
    left_join(dn |> mutate(cat = as.character(cat)), by = "cat") |>
    left_join(dt |> mutate(cat = as.character(cat)), by = "cat") |>
    left_join(cob |> mutate(cat = as.character(cat)), by = "cat") |>
    left_join(nn, by = "cat") |>
    mutate(bloco = bloco_de[[v]], variavel = rotulo_var[[v]], var_id = v,
           categoria   = cat,
           geral       = ic_fmt(geral_p, geral_lo, geral_hi),
           nao_testada = ic_fmt(nt_p, nt_lo, nt_hi),
           testada     = ic_fmt(t_p, t_lo, t_hi),
           cobertura   = ic_fmt(cob_p, cob_lo, cob_hi),
           p_valor     = ifelse(row_number() == 1, format.pval(pval, digits = 2, eps = 1e-4), "")) |>
    select(bloco, variavel, var_id, categoria, n, geral, nao_testada, testada, cobertura, p_valor)
}

tab1_cat <- bind_rows(lapply(vars_tab, linha_cat))

# idade contínua: média (IC) geral e por desfecho + p (teste t do desenho)
m_g <- sd |> summarise(m = survey_mean(idade, vartype = "ci"))
m_d <- sd |> group_by(desfecho_f) |> summarise(m = survey_mean(idade, vartype = "ci"))
p_idade <- svyttest(idade ~ desfecho_f, des)$p.value
linha_idade <- tibble(
  bloco = "Sociodemográfico", variavel = "Idade (anos, média)", var_id = "idade",
  categoria = "média (IC95%)", n = sum(!is.na(dc$idade)),
  geral = sprintf("%.1f (%.1f–%.1f)", m_g$m, m_g$m_low, m_g$m_upp),
  nao_testada = sprintf("%.1f (%.1f–%.1f)", m_d$m[m_d$desfecho_f == "Não testada"],
                        m_d$m_low[m_d$desfecho_f == "Não testada"], m_d$m_upp[m_d$desfecho_f == "Não testada"]),
  testada = sprintf("%.1f (%.1f–%.1f)", m_d$m[m_d$desfecho_f == "Testada"],
                    m_d$m_low[m_d$desfecho_f == "Testada"], m_d$m_upp[m_d$desfecho_f == "Testada"]),
  cobertura = "", p_valor = format.pval(p_idade, digits = 2, eps = 1e-4))

ordem_bloco <- c("Sociodemográfico", "Informação", "Reprodutivo", "Acesso", "Pós-exposição (serviço)")
tabela1 <- bind_rows(linha_idade, tab1_cat) |>
  mutate(bloco = factor(bloco, levels = ordem_bloco)) |>
  arrange(bloco, factor(var_id, levels = c("idade", vars_tab))) |>
  select(-var_id)
write_csv(tabela1, file.path(path_tabs, "03_tabela1_desfecho.csv"))

cat("Etapa 03 concluída.\n")
cat(sprintf("  Amostra principal: %d mulheres, %d países, %d PSU.\n",
            nrow(dc), n_distinct(dc$pais), n_distinct(dc$psu_id)))
cat(sprintf("  Prevalência pooled (w_igual): %.1f%% (IC95%% %.1f–%.1f).\n",
            pct1(prev_global$p), pct1(prev_global$p_low), pct1(prev_global$p_upp)))
cat("  Tabelas: 03_resumo_amostra.csv, 03_desfecho_por_pais.csv, 03_tabela1_desfecho.csv\n")
