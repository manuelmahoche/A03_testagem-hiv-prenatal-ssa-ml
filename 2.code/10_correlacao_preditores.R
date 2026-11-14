# ============================================================================
# P06 - Testagem de HIV no pré-natal na África Subsaariana
# Etapa 10 - Matriz de correlação/colinearidade dos preditores (análogo ao ESM8 da referência)
# Autor: Manuel Mahoche
# Entrada: 1.data/process/ds_clean.parquet
# Saída:   3.output/figs/10_fig_correlacao.png     mapa de calor das correlações (one-hot)
#          3.output/tabs/10_pares_colineares.csv   pares com |r| elevado (diagnóstico)
#
# Correlação de Pearson entre as variáveis indicadoras (one-hot) dos preditores, no
# Estrato A. Serve de diagnóstico de colinearidade (explica, p.ex., o rank-deficiency da
# regressão logística: os níveis "Sem parceiro" de educ/ocup do parceiro e decisão de saúde
# são redundantes com o estado conjugal "Fora de união").
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(forcats); library(here)
  library(ggplot2)
})
source(here("2.code", "R", "config.R"))
source(here("2.code", "R", "funcoes.R"))
dir.create(path_figs, showWarnings = FALSE, recursive = TRUE)

pred <- setdiff(vars_preditores, "regiao_id")
dc <- ler_parquet(file.path(path_process, "ds_clean.parquet")) |>
  filter(amostra == "principal", estrato == "A") |>
  select(all_of(pred))
for (v in pred) { x <- as.factor(dc[[v]]); if (anyNA(x)) { x <- addNA(x); levels(x)[is.na(levels(x))] <- "(sem dado)" }; dc[[v]] <- x }

X <- model.matrix(~ . - 1, data = dc)
# remove colunas constantes (sem variação) para a correlação
X <- X[, apply(X, 2, function(z) stats::sd(z) > 0), drop = FALSE]
R <- stats::cor(X)

# Pares com |r| elevado (exclui diagonal), para diagnóstico de colinearidade
rl <- as.data.frame(as.table(R)) |>
  filter(as.character(Var1) < as.character(Var2), abs(Freq) >= 0.4) |>
  arrange(desc(abs(Freq))) |>
  rename(variavel_1 = Var1, variavel_2 = Var2, r = Freq) |>
  mutate(r = round(r, 3))
write_csv(rl, file.path(path_tabs, "10_pares_colineares.csv"))

# Mapa de calor (paleta divergente: vermelho negativo, cinza no zero, azul positivo)
ord <- colnames(X)
long <- as.data.frame(as.table(R)) |>
  mutate(Var1 = factor(Var1, levels = ord), Var2 = factor(Var2, levels = rev(ord)))
p <- ggplot(long, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient2(low = "#B2182B", mid = "grey95", high = "#2166AC",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(title = "Correlação entre indicadores dos preditores (one-hot, Estrato A)",
       subtitle = "Diagnóstico de colinearidade; azul = positiva, vermelho = negativa",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
        axis.text.y = element_text(size = 5.5),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(color = "grey35", size = 9),
        legend.position = "right")
ggsave(file.path(path_figs, "10_fig_correlacao.png"), p, width = 9, height = 8, dpi = 300)

cat("Etapa 10 concluída.\n")
cat(sprintf("  %d colunas one-hot; %d pares com |r| >= 0.4.\n", ncol(X), nrow(rl)))
cat("  Top pares colineares:\n"); print(utils::head(rl, 10))
cat("  Figura: 10_fig_correlacao.png; Tabela: 10_pares_colineares.csv\n")
