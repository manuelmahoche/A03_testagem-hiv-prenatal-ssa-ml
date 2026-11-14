# ============================================================================
# P06 - Testagem de HIV no pré-natal na África Subsaariana
# Etapa 01 - População elegível, desfecho e controlo de qualidade (STATcompiler)
# Autor: Manuel Mahoche
# Entrada: 1.data/raw/dhs_ir_anc.parquet
# Saída:   1.data/process/ds_raw.parquet        (mulheres elegíveis + desfecho)
#          1.data/process/country_meta.parquet  (metadados e estrato por país)
#          3.output/tabs/01_qc_desfecho_statcompiler.csv
#          3.output/tabs/01_fluxo_amostra.csv
#
# Definições (DHS Guide to Statistics, indicador HA_CATH_W_ATR):
#   População: mulheres 15-49 anos com nascido vivo nos 24 meses anteriores à entrevista
#              (idade do filho mais recente = b19_01, ou v008 - b3_01 nas fases sem b19),
#              a quem foi aplicado o módulo de HIV (v750 não ausente; resolve as subamostras
#              do Gabão, Benim e Quénia) e fora dos campos de refugiados (Malawi 2024, scamps = 1).
#   Desfecho:  testada para HIV e recebeu o resultado no pré-natal (v841 = 1) ou entre a ida
#              para o parto e o nascimento (v841a = 1). Mulheres sem pré-natal entram com 0.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(here)
})
source(here("2.code", "R", "config.R"))
source(here("2.code", "R", "funcoes.R"))

ir <- ler_parquet(file.path(path_raw, "dhs_ir_anc.parquet"))

# ----------------------------- Elegibilidade ---------------------------------
ir <- ir |>
  mutate(
    idade_filho_m = coalesce(as.numeric(b19_01), as.numeric(v008) - as.numeric(b3_01)),
    e1_nasc24m    = !is.na(idade_filho_m) & idade_filho_m < 24,
    e2_modulo_hiv = !is.na(v750),
    e3_fora_campo = is.na(scamps) | scamps != 1
  )

fluxo <- ir |>
  group_by(pais) |>
  summarise(
    n0_todas        = n(),
    n1_nasc_24m     = sum(e1_nasc24m),
    n2_modulo_hiv   = sum(e1_nasc24m & e2_modulo_hiv),
    n3_fora_campos  = sum(e1_nasc24m & e2_modulo_hiv & e3_fora_campo),
    .groups = "drop"
  )

ds <- ir |>
  filter(e1_nasc24m, e2_modulo_hiv, e3_fora_campo) |>
  select(-starts_with("e1_"), -starts_with("e2_"), -starts_with("e3_"))

# ----------------------------- Desfecho --------------------------------------
ds <- ds |>
  mutate(
    hiv_anc_resultado = as.integer(v841 %in% 1 | v841a %in% 1),   # desfecho principal
    hiv_anc_testada   = as.integer(v840 %in% 1 | v840a %in% 1),   # testada (com ou sem resultado)
    hiv_teste_no_anc  = as.integer(v841 %in% 1),                  # resultado no pré-natal
    hiv_teste_no_parto = as.integer(!(v841 %in% 1) & v841a %in% 1),
    peso = as.numeric(v005) / 1e6
  )

# ----------------------------- Controlo de qualidade -------------------------
qc <- ds |>
  group_by(pais) |>
  summarise(
    n_elegiveis = n(),
    prev_pct    = round(100 * prop_ponderada(hiv_anc_resultado, peso), 1),
    testada_pct = round(100 * prop_ponderada(hiv_anc_testada, peso), 1),
    .groups = "drop"
  ) |>
  left_join(statcompiler, by = "pais") |>
  mutate(dif_pp = round(prev_pct - stat_pct, 1),
         dif_n  = n_elegiveis - stat_n,
         qc_ok  = abs(dif_pp) <= tolerancia_qc_pp)

falhas <- qc |> filter(!qc_ok)
verificar(nrow(falhas) == 0,
          paste0("Desfecho fora da tolerância do STATcompiler em: ",
                 paste(falhas$pais, collapse = ", "),
                 ". Rever o mapeamento v840-v841a antes de continuar."))

# ----------------------------- Metadados por país ----------------------------
meta <- inqueritos |>
  left_join(
    ds |> group_by(pais) |>
      summarise(fase = first(fase), ano_min = min(v007), ano_max = max(v007),
                n_psu = n_distinct(v021), .groups = "drop"),
    by = "pais") |>
  left_join(select(qc, pais, n_elegiveis, prev_pct, stat_pct, dif_pp), by = "pais") |>
  mutate(estrato = if_else(prev_pct < corte_estrato_pct, "A", "B"))

ds <- ds |> left_join(select(meta, pais, nome_pais, regiao_ssa, amostra, estrato), by = "pais")

# ----------------------------- Gravação --------------------------------------
gravar_parquet(ds,   file.path(path_process, "ds_raw.parquet"))
gravar_parquet(meta, file.path(path_process, "country_meta.parquet"))
write_csv(qc,    file.path(path_tabs, "01_qc_desfecho_statcompiler.csv"))
write_csv(fluxo, file.path(path_tabs, "01_fluxo_amostra.csv"))

cat("Etapa 01 concluída:", nrow(ds), "mulheres elegíveis;",
    sum(qc$qc_ok), "de", nrow(qc), "países dentro da tolerância do STATcompiler\n")
print(as.data.frame(select(qc, pais, n_elegiveis, prev_pct, stat_pct, dif_pp)))
