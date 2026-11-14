# ============================================================================
# P06 - Etapa 15 - Sensibilidade (janela 2015-17) e recalibração no Estrato B
# Autor: Manuel Mahoche
# Entrada: ds_clean.parquet, particoes.parquet, 06_modelos_finais.rds
# Saída:   3.output/tabs/15_sensibilidade_2015_17.csv   AUC nos 5 países da janela ampliada
#          3.output/tabs/15_recalibracao_B.csv          calibração no Estrato B, antes/depois
#
# (1) Sensibilidade temporal: o modelo de transportabilidade (sem país), treinado em TODO o
#     Estrato A (janela >=2018), é aplicado aos inquéritos de 2015-17 reservados (BJ, BU, UG,
#     ZA, ZW). AUC por país testa a generalização para uma janela temporal mais antiga.
# (2) Recalibração: no Estrato B, comparam-se calibração/Brier do modelo bruto e após
#     recalibração logística (Platt), mantendo a discriminação (AUC) inalterada.
# ============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(here); library(xgboost); library(pROC) })
source(here("2.code", "R", "config.R")); source(here("2.code", "R", "funcoes.R")); set.seed(semente)

pred <- setdiff(vars_preditores, "regiao_id")
dc <- ler_parquet(file.path(path_process, "ds_clean.parquet")) |>
  select(id_mulher, pais, nome_pais, amostra, estrato, all_of(pred), y = hiv_anc_resultado)
for (v in pred) { x <- as.factor(dc[[v]]); if (anyNA(x)) { x <- addNA(x); levels(x)[is.na(levels(x))] <- "(sem dado)" }; dc[[v]] <- x }
X <- model.matrix(~ . - 1, data = dc[pred]); y <- dc$y
idx_trainA <- which(dc$amostra == "principal" & dc$estrato == "A")

hp <- readRDS(file.path("3.output", "mod", "06_modelos_finais.rds"))$hp_xgb
m <- xgb.train(list(objective = "binary:logistic", eval_metric = "auc", max_depth = hp$max_depth,
                    eta = hp$eta, subsample = 0.8, colsample_bytree = 0.8, nthread = 0),
               xgb.DMatrix(X[idx_trainA, ], label = y[idx_trainA]), nrounds = hp$nrounds, verbose = 0)
dc$p <- as.numeric(predict(m, X))
auc1 <- function(yy, pp) tryCatch(round(as.numeric(pROC::auc(pROC::roc(yy, pp, quiet = TRUE, direction = "<"))), 3), error = function(e) NA_real_)

# (1) Sensibilidade 2015-17
sens <- dc |> filter(amostra == "sensibilidade") |> group_by(pais, nome_pais) |>
  summarise(n = n(), prev = round(100 * mean(y), 1), auc = auc1(y, p), .groups = "drop") |> arrange(desc(auc))
sens <- bind_rows(sens, dc |> filter(amostra == "sensibilidade") |>
                    summarise(pais = "MÉDIA", nome_pais = "Média (5 países 2015-17)", n = n(),
                              prev = round(100 * mean(y), 1), auc = auc1(y, p)))
write_csv(sens, file.path(path_tabs, "15_sensibilidade_2015_17.csv"))

# (2) Recalibração no Estrato B (Platt: y ~ logit(p))
B <- dc |> filter(amostra == "principal", estrato == "B")
clamp <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)
cal_fit <- glm(y ~ qlogis(clamp(p)), family = binomial, data = B)
B$p_recal <- as.numeric(predict(cal_fit, type = "response"))
cal_metrics <- function(yy, pp) {
  pp <- clamp(pp)
  int <- coef(glm(yy ~ offset(qlogis(pp)), family = binomial))[1]
  sl  <- coef(glm(yy ~ qlogis(pp), family = binomial))[2]
  c(auc = auc1(yy, pp), brier = round(mean((pp - yy)^2), 3),
    cal_int = round(unname(int), 2), cal_slope = round(unname(sl), 2))
}
recal <- bind_rows(
  tibble(versao = "Bruto (modelo do Estrato A)", !!!cal_metrics(B$y, B$p)),
  tibble(versao = "Recalibrado (Platt no Estrato B)", !!!cal_metrics(B$y, B$p_recal)))
write_csv(recal, file.path(path_tabs, "15_recalibracao_B.csv"))

cat("Etapa 15 concluída.\n")
cat("  Sensibilidade 2015-17 (AUC por país):\n"); print(as.data.frame(sens))
cat("  Recalibração no Estrato B:\n"); print(as.data.frame(recal))
