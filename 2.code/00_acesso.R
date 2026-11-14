# ============================================================================
# P06 - Testagem de HIV no pré-natal na África Subsaariana
# Etapa 00 - Leitura dos IR (Parquet) e empilhamento das variáveis do estudo
# Autor: Manuel Mahoche
# Entrada: <path_ir>/XXIRnnDT.parquet (um por inquérito, ver config.R)
# Saída:   1.data/raw/dhs_ir_anc.parquet (todas as mulheres, variáveis selecionadas)
#          3.output/tabs/00_inventario_variaveis.csv (disponibilidade por inquérito)
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tidyr); library(readr); library(here)
})
source(here("2.code", "R", "config.R"))
source(here("2.code", "R", "funcoes.R"))

verificar(dir.exists(path_ir) || identical(getOption("p06.io"), "csv"),
          paste0("Pasta dos IR não encontrada: ", path_ir,
                 "\nDefina DHS_PARQUET_IR ou ajuste path_ir em 2.code/R/config.R."))

ler_inquerito <- function(pais, ficheiro, ...) {
  f <- file.path(path_ir, paste0(ficheiro, ".parquet"))
  d <- ler_parquet(f, cols = vars_ler)
  # remove rótulos haven (se existirem) para empilhar países sem conflitos de classe
  d <- mutate(d, across(everything(), ~ { attributes(.x) <- NULL; .x }))
  # texto: formato varia por país -> forçar a character antes de empilhar
  for (v in intersect(vars_texto, names(d))) d[[v]] <- as.character(d[[v]])
  falta <- setdiff(vars_ler, names(d))
  for (v in falta) d[[v]] <- if (v %in% vars_texto) NA_character_ else NA_real_
  d |>
    mutate(pais = pais, ficheiro = ficheiro, fase = substr(as.character(v000), 3, 3)) |>
    select(pais, ficheiro, fase, all_of(vars_ler))
}

ir <- pmap(inqueritos, ler_inquerito) |> list_rbind()

# Inventário: % não ausente de cada variável, por inquérito (todas as mulheres)
inventario <- ir |>
  group_by(pais) |>
  summarise(across(all_of(vars_ler), ~ round(100 * mean(!is.na(.x)), 1)), n = n(), .groups = "drop")

gravar_parquet(ir, file.path(path_raw, "dhs_ir_anc.parquet"))
dir.create(path_tabs, showWarnings = FALSE, recursive = TRUE)
write_csv(inventario, file.path(path_tabs, "00_inventario_variaveis.csv"))

cat("Etapa 00 concluída:", nrow(ir), "mulheres em", n_distinct(ir$pais), "inquéritos\n")
