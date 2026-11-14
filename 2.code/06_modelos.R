# ============================================================================
# P06 - Testagem de HIV no pré-natal na África Subsaariana
# Etapa 06 - Modelos de aprendizado de máquina e transportabilidade
# Autor: Manuel Mahoche
# Entrada: 1.data/process/ds_clean.parquet, 1.data/process/particoes.parquet
# Saída:   3.output/tabs/06_hiperparametros.csv      hiperparâmetros afinados por CV
#          3.output/tabs/06_desempenho_interno.csv   métricas no teste interno (Estrato A)
#          3.output/tabs/06_desempenho_loco.csv      AUC por país retido (leave-one-country-out)
#          3.output/tabs/06_desempenho_externoB.csv  AUC em cada país do Estrato B
#          1.data/process/pred_teste_interno.parquet predições no teste interno (para Etapa 07)
#          3.output/mod/06_modelos_finais.rds         modelos finais (treinados no treino A)
#
# Desenho (ver Etapa 05):
#   - Modelagem no Estrato A; treino/teste bloqueados por PSU; 5 folds de CV (por PSU) p/ afinar.
#   - Preditores: vars_preditores SEM regiao_id (evita identificar o país e viciar o LOCO).
#     Pós-exposição fica de fora (mediadores). Ausências entram como nível "(sem dado)".
#   - Algoritmos: regressão logística, logística penalizada (glmnet), random forest (ranger),
#     gradient boosting (xgboost). Treino NÃO ponderado (prática padrão de predição).
#   - Afinação por CV (folds por PSU) só no treino; os hiperparâmetros escolhidos são reusados
#     nos reajustes do LOCO e do Estrato B (sem reafinar por fold).
#   - Transportabilidade: LOCO (16 países A) e teste externo no Estrato B (7 países).
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(here)
  library(glmnet); library(ranger); library(xgboost); library(pROC)
})
source(here("2.code", "R", "config.R"))
source(here("2.code", "R", "funcoes.R"))
set.seed(semente)

pred <- setdiff(vars_preditores, "regiao_id")

# ---------------------------------------------------------------------------
# Dados: junta ds_clean (desfecho + covariáveis) às partições, só amostra principal
# ---------------------------------------------------------------------------
dc <- ler_parquet(file.path(path_process, "ds_clean.parquet")) |>
  filter(amostra == "principal") |>
  select(id_mulher, pais, nome_pais, estrato, all_of(pred), y = hiv_anc_resultado)
pa <- ler_parquet(file.path(path_process, "particoes.parquet")) |>
  select(id_mulher, particao, cv_fold)
dat <- inner_join(dc, pa, by = "id_mulher")

# Ausências -> nível explícito "(sem dado)"; garante factor em todos os preditores
for (v in pred) {
  x <- as.factor(dat[[v]])
  if (anyNA(x)) { x <- addNA(x); levels(x)[is.na(levels(x))] <- "(sem dado)" }
  dat[[v]] <- x
}
dat$yf <- factor(dat$y, levels = c(0, 1))

# Matriz one-hot (para glmnet/xgboost), construída sobre TODAS as linhas -> colunas idênticas
X <- model.matrix(~ . - 1, data = dat[pred])
y <- dat$y

idx_tr   <- which(dat$particao == "treino")
idx_te   <- which(dat$particao == "teste")
idx_B    <- which(dat$estrato == "B")
idx_allA <- which(dat$estrato == "A")                      # treino+teste A (para o teste externo B)
fold_tr  <- dat$cv_fold[idx_tr]

cat(sprintf("Dados: treino=%d, teste interno=%d, Estrato B=%d; preditores=%d, colunas one-hot=%d\n",
            length(idx_tr), length(idx_te), length(idx_B), length(pred), ncol(X)))

# ---------------------------------------------------------------------------
# Métricas
# ---------------------------------------------------------------------------
auprc <- function(yv, p) {
  o <- order(p, decreasing = TRUE); yv <- yv[o]
  tp <- cumsum(yv); fp <- cumsum(1 - yv)
  prec <- tp / (tp + fp); rec <- tp / sum(yv)
  sum(prec * diff(c(0, rec)))
}
calib <- function(yv, p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6); lp <- qlogis(p)
  int <- tryCatch(unname(coef(glm(yv ~ offset(lp), family = binomial))[1]), error = function(e) NA_real_)
  sl  <- tryCatch(unname(coef(glm(yv ~ lp, family = binomial))[2]), error = function(e) NA_real_)
  c(cal_int = int, cal_slope = sl)
}
metrics <- function(yv, p) {
  a <- tryCatch(as.numeric(pROC::auc(pROC::roc(yv, p, quiet = TRUE, direction = "<"))),
                error = function(e) NA_real_)
  c(auc = a, auprc = auprc(yv, p), brier = mean((p - yv)^2), calib(yv, p))
}
auc_only <- function(yv, p) tryCatch(
  as.numeric(pROC::auc(pROC::roc(yv, p, quiet = TRUE, direction = "<"))), error = function(e) NA_real_)

# ---------------------------------------------------------------------------
# Definição dos algoritmos: fit(rows, hp) e predict em linhas novas
#   Cada learner devolve uma função que, dado o conjunto de treino (linhas) e os
#   hiperparâmetros, produz um preditor de probabilidades para novas linhas.
# ---------------------------------------------------------------------------
fit_lr <- function(rows, hp = NULL) {
  m <- glm(y ~ ., data = data.frame(y = dat$y[rows], X[rows, , drop = FALSE]), family = binomial)
  function(nr) as.numeric(predict(m, newdata = data.frame(X[nr, , drop = FALSE]), type = "response"))
}
fit_glmnet <- function(rows, hp) {
  m <- glmnet(X[rows, , drop = FALSE], y[rows], family = "binomial", alpha = hp$alpha, lambda = hp$lambda)
  function(nr) as.numeric(predict(m, X[nr, , drop = FALSE], type = "response"))
}
fit_rf <- function(rows, hp) {
  m <- ranger(x = dat[rows, pred, drop = FALSE], y = dat$yf[rows], probability = TRUE,
              num.trees = hp$num.trees, mtry = hp$mtry, min.node.size = hp$min.node.size,
              num.threads = 0, seed = semente)
  function(nr) predict(m, dat[nr, pred, drop = FALSE])$predictions[, "1"]
}
fit_xgb <- function(rows, hp) {
  dtr <- xgb.DMatrix(X[rows, , drop = FALSE], label = y[rows])
  m <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc",
                               max_depth = hp$max_depth, eta = hp$eta,
                               subsample = 0.8, colsample_bytree = 0.8, nthread = 0),
                 data = dtr, nrounds = hp$nrounds, verbose = 0)
  function(nr) as.numeric(predict(m, X[nr, , drop = FALSE]))
}

# ---------------------------------------------------------------------------
# Afinação por CV (folds por PSU) no conjunto de treino
# ---------------------------------------------------------------------------
cat("Afinando hiperparâmetros por CV (folds por PSU)...\n")
foldid <- fold_tr

# glmnet: alpha em {0.5, 1}; lambda por cv.glmnet (AUC)
tune_glmnet <- function() {
  best <- NULL
  for (al in c(0.5, 1)) {
    cv <- cv.glmnet(X[idx_tr, , drop = FALSE], y[idx_tr], family = "binomial",
                    foldid = foldid, type.measure = "auc", alpha = al)
    auc_cv <- max(cv$cvm)
    if (is.null(best) || auc_cv > best$auc) best <- list(alpha = al, lambda = cv$lambda.min, auc = auc_cv)
  }
  best
}

# CV-AUC genérico (out-of-fold) para ranger/xgboost dado um construtor
cv_auc_oof <- function(fit_fun, hp) {
  oof <- rep(NA_real_, length(idx_tr))
  for (k in sort(unique(foldid))) {
    val_local <- which(foldid == k)
    tr_rows   <- idx_tr[foldid != k]
    predictor <- fit_fun(tr_rows, hp)
    oof[val_local] <- predictor(idx_tr[val_local])
  }
  auc_only(y[idx_tr], oof)
}
tune_rf <- function() {
  p <- length(pred)
  grid <- expand.grid(mtry = unique(pmax(1, round(c(sqrt(p), p / 3)))),
                      min.node.size = c(10, 50), num.trees = 300)
  grid$auc <- NA_real_
  for (i in seq_len(nrow(grid))) grid$auc[i] <- cv_auc_oof(fit_rf, as.list(grid[i, ]))
  as.list(grid[which.max(grid$auc), ])
}
tune_xgb <- function() {
  # nrounds por xgb.cv (early stopping) num ponto central; grid em depth/eta
  folds <- split(seq_along(idx_tr), foldid)
  grid <- expand.grid(max_depth = c(3, 6), eta = c(0.1, 0.3))
  grid$nrounds <- NA_integer_; grid$auc <- NA_real_
  dall <- xgb.DMatrix(X[idx_tr, , drop = FALSE], label = y[idx_tr])
  for (i in seq_len(nrow(grid))) {
    cv <- xgb.cv(params = list(objective = "binary:logistic", eval_metric = "auc",
                               max_depth = grid$max_depth[i], eta = grid$eta[i],
                               subsample = 0.8, colsample_bytree = 0.8, nthread = 0),
                 data = dall, nrounds = 500, folds = folds, early_stopping_rounds = 25,
                 verbose = 0)
    # Deriva o melhor nº de iterações do log (robusto às variações de API do xgboost).
    log <- as.data.frame(cv$evaluation_log)
    auc_col <- grep("test.*auc.*mean", names(log), ignore.case = TRUE, value = TRUE)[1]
    best_it <- which.max(log[[auc_col]])
    grid$nrounds[i] <- best_it
    grid$auc[i] <- log[[auc_col]][best_it]
  }
  as.list(grid[which.max(grid$auc), ])
}

hp_glmnet <- tune_glmnet(); cat(sprintf("  glmnet: alpha=%.1f lambda=%.4g (CV-AUC %.3f)\n",
                                        hp_glmnet$alpha, hp_glmnet$lambda, hp_glmnet$auc))
hp_rf <- tune_rf(); cat(sprintf("  ranger: mtry=%d min.node.size=%d trees=%d (CV-AUC %.3f)\n",
                                hp_rf$mtry, hp_rf$min.node.size, hp_rf$num.trees, hp_rf$auc))
hp_xgb <- tune_xgb(); cat(sprintf("  xgboost: depth=%d eta=%.2f nrounds=%d (CV-AUC %.3f)\n",
                                  hp_xgb$max_depth, hp_xgb$eta, hp_xgb$nrounds, hp_xgb$auc))

hiper <- bind_rows(
  tibble(modelo = "glmnet",  hiperparametros = sprintf("alpha=%.1f; lambda=%.4g", hp_glmnet$alpha, hp_glmnet$lambda), cv_auc = round(hp_glmnet$auc, 3)),
  tibble(modelo = "ranger",  hiperparametros = sprintf("mtry=%d; min.node.size=%d; num.trees=%d", hp_rf$mtry, hp_rf$min.node.size, hp_rf$num.trees), cv_auc = round(hp_rf$auc, 3)),
  tibble(modelo = "xgboost", hiperparametros = sprintf("max_depth=%d; eta=%.2f; nrounds=%d", hp_xgb$max_depth, hp_xgb$eta, hp_xgb$nrounds), cv_auc = round(hp_xgb$auc, 3)))
write_csv(hiper, file.path(path_tabs, "06_hiperparametros.csv"))

# ---------------------------------------------------------------------------
# Modelos finais (treino A) e avaliação no teste interno
# ---------------------------------------------------------------------------
cat("Treinando modelos finais e avaliando no teste interno...\n")
learners <- list(
  lr      = function(rows) fit_lr(rows),
  glmnet  = function(rows) fit_glmnet(rows, hp_glmnet),
  ranger  = function(rows) fit_rf(rows, hp_rf),
  xgboost = function(rows) fit_xgb(rows, hp_xgb))

preditores_fin <- lapply(learners, function(L) L(idx_tr))
pred_te <- sapply(preditores_fin, function(pf) pf(idx_te))

interno <- lapply(names(preditores_fin), function(nm) {
  m <- metrics(y[idx_te], pred_te[, nm])
  tibble(modelo = nm, auc = round(m["auc"], 3), auprc = round(m["auprc"], 3),
         brier = round(m["brier"], 3), cal_int = round(m["cal_int"], 2), cal_slope = round(m["cal_slope"], 2))
}) |> bind_rows() |> arrange(desc(auc))
write_csv(interno, file.path(path_tabs, "06_desempenho_interno.csv"))

# Predições do teste interno (para a Etapa 07)
pred_tbl <- tibble(id_mulher = dat$id_mulher[idx_te], pais = dat$pais[idx_te],
                   nome_pais = dat$nome_pais[idx_te], desfecho = y[idx_te]) |>
  bind_cols(as_tibble(pred_te) |> rename_with(~ paste0("p_", .x)))
gravar_parquet(pred_tbl, file.path(path_process, "pred_teste_interno.parquet"))

dir.create(file.path("3.output", "mod"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(hp_glmnet = hp_glmnet, hp_rf = hp_rf, hp_xgb = hp_xgb,
             interno = interno), file.path("3.output", "mod", "06_modelos_finais.rds"))

# ---------------------------------------------------------------------------
# LOCO: treina em 15 países A, testa no país retido (AUC por modelo)
# ---------------------------------------------------------------------------
cat("LOCO (16 países)...\n")
paises_A <- sort(unique(dat$pais[idx_allA]))
loco <- lapply(paises_A, function(pc) {
  tr_rows <- idx_allA[dat$pais[idx_allA] != pc]
  te_rows <- idx_allA[dat$pais[idx_allA] == pc]
  row <- tibble(pais = pc, nome_pais = dat$nome_pais[te_rows][1], n_teste = length(te_rows))
  for (nm in names(learners)) {
    pf <- tryCatch(learners[[nm]](tr_rows), error = function(e) NULL)
    row[[paste0("auc_", nm)]] <- if (is.null(pf)) NA_real_ else round(auc_only(y[te_rows], pf(te_rows)), 3)
  }
  row
}) |> bind_rows()
loco <- bind_rows(loco, loco |> summarise(pais = "MÉDIA", nome_pais = "Média (16 países)", n_teste = sum(n_teste),
                                          across(starts_with("auc_"), ~ round(mean(.x, na.rm = TRUE), 3))))
write_csv(loco, file.path(path_tabs, "06_desempenho_loco.csv"))

# ---------------------------------------------------------------------------
# Teste externo no Estrato B: treina em TODO o Estrato A, prediz cada país B
# ---------------------------------------------------------------------------
cat("Teste externo (Estrato B)...\n")
pred_B <- lapply(learners, function(L) L(idx_allA))
paises_B <- sort(unique(dat$pais[idx_B]))
externoB <- lapply(paises_B, function(pc) {
  te_rows <- idx_B[dat$pais[idx_B] == pc]
  row <- tibble(pais = pc, nome_pais = dat$nome_pais[te_rows][1], n = length(te_rows),
                prev = round(100 * mean(y[te_rows]), 1))
  for (nm in names(pred_B)) row[[paste0("auc_", nm)]] <- round(auc_only(y[te_rows], pred_B[[nm]](te_rows)), 3)
  row
}) |> bind_rows()
externoB <- bind_rows(externoB, externoB |> summarise(pais = "MÉDIA", nome_pais = "Média (7 países B)",
                                                      n = sum(n), prev = round(mean(prev), 1),
                                                      across(starts_with("auc_"), ~ round(mean(.x, na.rm = TRUE), 3))))
write_csv(externoB, file.path(path_tabs, "06_desempenho_externoB.csv"))

cat("Etapa 06 concluída.\n")
cat("  Teste interno (AUC):\n"); print(as.data.frame(interno[, c("modelo", "auc", "auprc", "brier")]))
cat("  LOCO média e Estrato B média:\n")
print(as.data.frame(loco[loco$pais == "MÉDIA", grep("nome_pais|auc_", names(loco))]))
print(as.data.frame(externoB[externoB$pais == "MÉDIA", grep("nome_pais|auc_", names(externoB))]))
cat("  Tabelas em 3.output/tabs/06_*.csv; predições e modelos guardados.\n")
