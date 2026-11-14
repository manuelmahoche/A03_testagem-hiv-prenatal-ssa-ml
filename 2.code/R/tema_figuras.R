# ============================================================================
# P06 - Tema gráfico partilhado (qualidade de publicação) e paleta Okabe-Ito
# Uso: source(here("2.code","R","tema_figuras.R")) nos scripts de figuras
# ============================================================================
suppressPackageStartupMessages(library(ggplot2))

# Paleta qualitativa segura para daltonismo (Okabe-Ito)
cores_okabe <- c("#0072B2", "#D55E00", "#009E73", "#E69F00",
                 "#56B4E9", "#CC79A7", "#F0E442", "#000000")
cores_modelo <- c(xgboost = "#0072B2", ranger = "#D55E00", glmnet = "#009E73", lr = "#E69F00")

# Tema base de publicação
tema_pub <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      plot.title    = element_text(face = "bold", size = base + 1),
      plot.subtitle = element_text(color = "grey35", size = base - 2, margin = margin(b = 6)),
      plot.caption  = element_text(color = "grey45", size = base - 4, hjust = 0),
      axis.title    = element_text(size = base - 1, color = "grey20"),
      axis.text     = element_text(color = "grey30"),
      strip.text    = element_text(face = "bold", size = base - 2),
      legend.position = "top", legend.title = element_blank(),
      legend.text   = element_text(size = base - 2),
      plot.margin   = margin(12, 16, 10, 12))
}
