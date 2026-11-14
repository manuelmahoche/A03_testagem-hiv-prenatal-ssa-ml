# ============================================================================
# P06 - Testagem de HIV no pré-natal na África Subsaariana
# Etapa 07 - Calibração e interpretabilidade dos modelos
# Autor: Manuel Mahoche
# Entrada: ds_clean.parquet, particoes.parquet, pred_teste_interno.parquet,
#          3.output/mod/06_modelos_finais.rds (hiperparâmetros)
# Saída:   3.output/tabs/07_calibracao_resumo.csv   int/incl./Brier/ECE por modelo
#          3.output/tabs/07_calibracao_bins.csv      curva de calibração (decis) por modelo
#          3.output/tabs/07_importancia.csv          importância por variável (3 modelos)
#          3.output/figs/07_fig_calibracao.png       curvas de calibração
#          3.output/figs/07_fig_importancia_shap.png importância SHAP (xgboost)
#          3.output/figs/07_fig_transportabilidade.png AUC por país (LOCO + externo B)
#
# Calibração no teste interno (reaproveita pred_teste_interno.parquet da Etapa 06).
# Interpretabilidade: reajusta os modelos no treino com os hiperparâmetros afinados e
#   - xgboost: valores SHAP exatos via predict(predcontrib = TRUE) no teste;
#   - ranger:  importância por permutação (preditores em fator, já por variável);
#   - glmnet:  |coeficiente| agregado por variável.
# Paleta categórica: Okabe-Ito (segura para daltonismo).
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(forcats); library(here)
  library(glmnet); library(ranger); library(xgboost); library(ggplot2); library(scales)
})
source(here("2.code", "R", "config.R"))
source(here("2.code", "R", "funcoes.R"))
set.seed(semente)
dir.create(path_figs, showWarnings = FALSE, recursive = TRUE)

# Paleta Okabe-Ito por modelo e tema de publicação partilhado ----------------
source(here("2.code", "R", "tema_figuras.R"))
rotulo_modelo <- c(lr = "Regressão logística", glmnet = "Logística penalizada (glmnet)",
                   ranger = "Random forest (ranger)", xgboost = "Gradient boosting (xgboost)")
tema_p06 <- tema_pub(12)

# ---------------------------------------------------------------------------
# 1) Calibração (teste interno)
# ---------------------------------------------------------------------------
pt <- ler_parquet(file.path(path_process, "pred_teste_interno.parquet"))
mods <- c("lr", "glmnet", "ranger", "xgboost")

long <- pt |>
  select(desfecho, starts_with("p_")) |>
  pivot_longer(starts_with("p_"), names_to = "modelo", values_to = "p") |>
  mutate(modelo = sub("^p_", "", modelo))

# Curva por decis de probabilidade prevista
bins <- long |> group_by(modelo) |>
  mutate(bin = ntile(p, 10)) |>
  group_by(modelo, bin) |>
  summarise(pred_media = mean(p), obs_frac = mean(desfecho), n = n(), .groups = "drop")
write_csv(bins, file.path(path_tabs, "07_calibracao_bins.csv"))

# Resumo: intercepto, inclinação, Brier, ECE
calib_resumo <- long |> group_by(modelo) |>
  summarise(
    brier = mean((p - desfecho)^2),
    cal_int = {pp <- pmin(pmax(p, 1e-6), 1 - 1e-6); coef(glm(desfecho ~ offset(qlogis(pp)), family = binomial))[1]},
    cal_slope = {pp <- pmin(pmax(p, 1e-6), 1 - 1e-6); coef(glm(desfecho ~ qlogis(pp), family = binomial))[2]},
    .groups = "drop") |>
  left_join(bins |> group_by(modelo) |> summarise(ece = weighted.mean(abs(obs_frac - pred_media), n), .groups = "drop"),
            by = "modelo") |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  arrange(brier)
write_csv(calib_resumo, file.path(path_tabs, "07_calibracao_resumo.csv"))

p_cal <- ggplot(bins, aes(pred_media, obs_frac, color = modelo)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  scale_color_manual(values = cores_modelo, labels = rotulo_modelo) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Calibração no teste interno (Estrato A)",
       subtitle = "Proporção observada de testagem vs. probabilidade prevista, por decis",
       x = "Probabilidade prevista", y = "Proporção observada") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  coord_equal() + tema_p06 + theme(legend.text = element_text(size = 9))
ggsave(file.path(path_figs, "07_fig_calibracao.png"), p_cal, width = 6.5, height = 6.5, dpi = 300)

# ---------------------------------------------------------------------------
# 2) Reconstrói dados/matriz (igual à Etapa 06) e reajusta modelos no treino
# ---------------------------------------------------------------------------
pred <- setdiff(vars_preditores, "regiao_id")
dc <- ler_parquet(file.path(path_process, "ds_clean.parquet")) |>
  filter(amostra == "principal") |>
  select(id_mulher, pais, nome_pais, estrato, all_of(pred), y = hiv_anc_resultado)
pa <- ler_parquet(file.path(path_process, "particoes.parquet")) |> select(id_mulher, particao, cv_fold)
dat <- inner_join(dc, pa, by = "id_mulher")
for (v in pred) { x <- as.factor(dat[[v]]); if (anyNA(x)) { x <- addNA(x); levels(x)[is.na(levels(x))] <- "(sem dado)" }; dat[[v]] <- x }
dat$yf <- factor(dat$y, levels = c(0, 1))
X <- model.matrix(~ . - 1, data = dat[pred]); y <- dat$y
idx_tr <- which(dat$particao == "treino"); idx_te <- which(dat$particao == "teste")

hp <- readRDS(file.path("3.output", "mod", "06_modelos_finais.rds"))

# Mapeia cada coluna one-hot à variável de origem (prefixo mais longo)
col2var <- vapply(colnames(X), function(cn) { h <- pred[startsWith(cn, pred)]; h[which.max(nchar(h))] }, character(1))

# -- xgboost: SHAP exato no teste --
m_xgb <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc",
                                 max_depth = hp$hp_xgb$max_depth, eta = hp$hp_xgb$eta,
                                 subsample = 0.8, colsample_bytree = 0.8, nthread = 0),
                   data = xgb.DMatrix(X[idx_tr, ], label = y[idx_tr]), nrounds = hp$hp_xgb$nrounds, verbose = 0)
shap <- predict(m_xgb, X[idx_te, ], predcontrib = TRUE)        # n x (p+1), última coluna = BIAS
shap <- shap[, colnames(shap) != "BIAS", drop = FALSE]
shap_col <- colMeans(abs(shap))                                 # média |SHAP| por coluna one-hot
shap_var <- tapply(shap_col, col2var[colnames(shap)], sum)      # agrega por variável

# -- ranger: importância por permutação (por variável) --
m_rf <- ranger(x = dat[idx_tr, pred], y = dat$yf[idx_tr], probability = TRUE,
               num.trees = hp$hp_rf$num.trees, mtry = hp$hp_rf$mtry, min.node.size = hp$hp_rf$min.node.size,
               importance = "permutation", num.threads = 0, seed = semente)
imp_rf <- ranger::importance(m_rf)

# -- glmnet: |coeficiente| agregado por variável --
m_gl <- glmnet(X[idx_tr, ], y[idx_tr], family = "binomial", alpha = hp$hp_glmnet$alpha, lambda = hp$hp_glmnet$lambda)
beta <- as.matrix(coef(m_gl)); beta <- beta[rownames(beta) != "(Intercept)", 1]
beta_var <- tapply(abs(beta), col2var[names(beta)], sum)

# as.numeric() tira os atributos de array que vêm do tapply (senão viram coluna-matriz no tibble)
norm100 <- function(v) { v <- as.numeric(v); round(100 * v / sum(v, na.rm = TRUE), 1) }
importancia <- tibble(variavel = pred) |>
  mutate(shap_xgb = norm100(shap_var[variavel]),
         imp_ranger = norm100(imp_rf[variavel]),
         absbeta_glmnet = norm100(beta_var[variavel]),
         media_rel = round(rowMeans(cbind(shap_xgb, imp_ranger, absbeta_glmnet), na.rm = TRUE), 1)) |>
  arrange(desc(shap_xgb))
write_csv(importancia, file.path(path_tabs, "07_importancia.csv"))

# Figura: top preditores por SHAP (xgboost)
rot_pred <- c(idade_grp = "Faixa etária", uniao = "Estado conjugal", residencia = "Residência",
              chefe_sexo = "Sexo do chefe", riqueza = "Riqueza", escolaridade = "Escolaridade",
              ocupacao = "Ocupação", educ_parceiro = "Escol. parceiro", ocup_parceiro = "Ocup. parceiro",
              media_semanal = "Exposição a média", internet_12m = "Usou internet", conhece_ptv = "Conhece PTV",
              estigma = "Estigma", paridade = "Paridade", idade_1parto = "Idade 1º parto",
              intencao_gravidez = "Intenção gravidez", desejo_filhos = "Desejo filhos",
              distancia_problema = "Distância problema", seguro_saude = "Seguro de saúde")
imp_fig <- importancia |> mutate(lbl = fct_reorder(rot_pred[variavel], shap_xgb))
p_imp <- ggplot(imp_fig, aes(shap_xgb, lbl)) +
  geom_col(fill = cores_modelo["xgboost"], width = 0.72) +
  geom_text(aes(label = sprintf("%.1f", shap_xgb)), hjust = -0.2, size = 3, color = "grey25") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Importância dos preditores (SHAP, xgboost)",
       subtitle = "Contribuição relativa (%) para a predição no teste interno",
       x = "Importância relativa (%)", y = NULL) + tema_p06 +
  theme(legend.position = "none", panel.grid.major.y = element_blank())
ggsave(file.path(path_figs, "07_fig_importancia_shap.png"), p_imp, width = 7, height = 5.5, dpi = 300)

# ---------------------------------------------------------------------------
# 3) Figura de transportabilidade: AUC por país (LOCO em A, externo em B)
# ---------------------------------------------------------------------------
loco <- read_csv(file.path(path_tabs, "06_desempenho_loco.csv"), show_col_types = FALSE) |> filter(pais != "MÉDIA")
extB <- read_csv(file.path(path_tabs, "06_desempenho_externoB.csv"), show_col_types = FALSE) |> filter(pais != "MÉDIA")
auc_interno <- read_csv(file.path(path_tabs, "06_desempenho_interno.csv"), show_col_types = FALSE) |>
  filter(modelo == "xgboost") |> pull(auc)

transp <- bind_rows(
  loco |> transmute(nome_pais, auc = auc_xgboost, tipo = "Estrato A (LOCO)"),
  extB |> transmute(nome_pais, auc = auc_xgboost, tipo = "Estrato B (externo)")) |>
  mutate(nome_pais = fct_reorder(nome_pais, auc))
p_tr <- ggplot(transp, aes(auc, nome_pais, color = tipo)) +
  geom_vline(xintercept = auc_interno, linetype = "dashed", color = "grey55") +
  geom_point(size = 2.6) +
  scale_color_manual(values = c("Estrato A (LOCO)" = "#0072B2", "Estrato B (externo)" = "#D55E00")) +
  scale_x_continuous(limits = c(0.5, 0.92), breaks = seq(0.5, 0.9, 0.1)) +
  labs(title = "Transportabilidade entre países (xgboost)",
       subtitle = sprintf("AUC no país retido (LOCO) e no Estrato B (externo); tracejado = teste interno (%.2f)", auc_interno),
       x = "AUC", y = NULL) + tema_p06
ggsave(file.path(path_figs, "07_fig_transportabilidade.png"), p_tr, width = 7, height = 6.5, dpi = 300)

cat("Etapa 07 concluída.\n")
cat("  Calibração (resumo):\n"); print(as.data.frame(calib_resumo))
cat("  Top 8 preditores (SHAP xgboost):\n"); print(as.data.frame(head(importancia, 8)))
cat("  Figuras: 07_fig_calibracao.png, 07_fig_importancia_shap.png, 07_fig_transportabilidade.png\n")
