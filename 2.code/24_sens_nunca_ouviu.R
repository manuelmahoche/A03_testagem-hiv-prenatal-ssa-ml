# ============================================================================
# P06 - Etapa 24 - Sensibilidade: mulheres que nunca ouviram falar de HIV (zero estrutural do desfecho)
# Autor: Manuel Mahoche
# Motivo: no DHS, quem nunca ouviu falar de AIDS (v751 = 0) salta as perguntas sobre o teste (v841 ausente) e
# conta como não testada (cobertura 0% por construção). A categoria "Não ouviu falar de HIV" do estigma é por
# isso um preditor determinístico. Aqui repete-se a modelagem SÓ entre quem ouviu falar de HIV.
# Saída: tabs/24_sens_nunca_ouviu.csv (AUC interno, LOCO, Estrato B, 2015-17; modelo completo vs restrito)
#        tabs/24_importancia_restrita.csv (SHAP por variável, modelo restrito, teste interno)
# ============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr); library(here); library(xgboost) })
source(here("2.code", "R", "config.R")); source(here("2.code", "R", "funcoes.R")); set.seed(semente)
pred <- setdiff(vars_preditores, "regiao_id")
dc <- ler_parquet(file.path(path_process, "ds_clean.parquet")) |>
  select(id_mulher, pais, psu_id, amostra, estrato, all_of(pred), y = hiv_anc_resultado)
pa <- ler_parquet(file.path(path_process, "particoes.parquet")) |> select(id_mulher, particao)
dc <- left_join(dc, pa, by = "id_mulher")
for (v in pred) { x <- as.factor(dc[[v]]); if (anyNA(x)) { x <- addNA(x); levels(x)[is.na(levels(x))] <- "(sem dado)" }; dc[[v]] <- x }
heard <- dc$estigma != "Não ouviu falar de HIV"; y <- dc$y
X <- model.matrix(~ . - 1, data = dc[pred]); X <- X[, colSums(X[dc$amostra == "principal", , drop = FALSE]) > 0, drop = FALSE]
hp <- readRDS(file.path("3.output", "mod", "06_modelos_finais.rds"))$hp_xgb
auc <- function(yy, p) { n1 <- sum(yy == 1); n0 <- length(yy) - n1; if (n1 == 0 || n0 == 0) NA else (sum(rank(p)[yy == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
fit <- function(rows) { m <- xgb.train(list(objective = "binary:logistic", eval_metric = "auc", max_depth = hp$max_depth, eta = hp$eta,
                        subsample = 0.8, colsample_bytree = 0.8, nthread = 0), xgb.DMatrix(X[rows, , drop = FALSE], label = y[rows]), nrounds = hp$nrounds, verbose = 0); m }
pr <- function(m, rows) as.numeric(predict(m, X[rows, , drop = FALSE]))
A <- which(dc$amostra == "principal" & dc$estrato == "A"); B <- which(dc$amostra == "principal" & dc$estrato == "B"); S <- which(dc$amostra == "sensibilidade")
tr <- A[dc$particao[A] == "treino"]; te <- A[dc$particao[A] == "teste"]
medias <- function(rows, p) mean(sapply(split(seq_along(rows), dc$pais[rows]), function(i) auc(y[rows][i], p[i])))
res <- list()
for (modo in c("Modelo completo (treino com todas)", "Modelo restrito a quem ouviu falar de HIV")) {
  restr <- modo != "Modelo completo (treino com todas)"
  keep <- function(r) if (restr) r[heard[r]] else r
  m_in <- fit(keep(tr))
  teh <- te[heard[te]]
  loco <- sapply(sort(unique(dc$pais[A])), function(pc) { trn <- keep(A[dc$pais[A] != pc]); tst <- A[dc$pais[A] == pc & heard[A]]; auc(y[tst], pr(fit(trn), tst)) })
  mA <- fit(keep(A)); Bh <- B[heard[B]]; Sh <- S[heard[S]]
  res[[modo]] <- tibble(modelo = modo, avaliado_em = "so quem ouviu falar de HIV",
    auc_interno = auc(y[teh], pr(m_in, teh)), auc_loco_media = mean(loco), auc_B_media = medias(Bh, pr(mA, Bh)),
    auc_2015_17_media = medias(Sh, pr(mA, Sh)), n_interno = length(teh))
  if (restr) {
    contrib <- predict(m_in, X[teh, , drop = FALSE], predcontrib = TRUE); contrib <- contrib[, colnames(contrib) %in% colnames(X), drop = FALSE]
    var_de <- vapply(colnames(contrib), function(cn) { v <- pred[startsWith(cn, pred)]; v[which.max(nchar(v))] }, character(1))
    imp <- tapply(colMeans(abs(contrib)), var_de, sum); imp <- round(100 * imp / sum(imp), 1)
    write_csv(tibble(variavel = names(imp), shap_pct = as.numeric(imp)) |> arrange(desc(shap_pct)), file.path(path_tabs, "24_importancia_restrita.csv"))
  }
}
out <- bind_rows(res) |> mutate(across(where(is.numeric) & !n_interno, ~ round(.x, 3)))
write_csv(out, file.path(path_tabs, "24_sens_nunca_ouviu.csv")); print(as.data.frame(out))
print(read_csv(file.path(path_tabs, "24_importancia_restrita.csv"), show_col_types = FALSE) |> head(6))
