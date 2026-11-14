# ============================================================================
# P06 - Corre a preparação dos dados (Etapas 00 a 02) por ordem
# ============================================================================
etapas <- c("00_acesso.R", "01_elegibilidade_desfecho.R", "02_limpeza.R")
for (e in etapas) {
  cat("\n==========", e, "==========\n")
  source(here::here("2.code", e), local = new.env())
}
cat("\nPreparação concluída. Ver 3.output/tabs/01_qc_desfecho_statcompiler.csv\n")
