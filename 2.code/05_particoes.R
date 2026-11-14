# ============================================================================
# P06 - Testagem de HIV no pré-natal na África Subsaariana
# Etapa 05 - Partições de validação para a modelagem (split por PSU, CV e LOCO)
# Autor: Manuel Mahoche
# Entrada: 1.data/process/ds_clean.parquet
# Saída:   1.data/process/particoes.parquet        uma linha por mulher (amostra principal)
#          3.output/tabs/05_particoes_resumo.csv    balanço por país (n e prevalência)
#          3.output/tabs/05_loco_folds.csv          tamanho de cada fold LOCO (16 países A)
#
# Desenho de validação (evita fuga por conglomerado — mulheres do mesmo PSU nunca
# ficam em lados diferentes):
#   - Modelagem no ESTRATO A (16 países, cobertura <80%).
#   - Holdout: PSU aleatórios para treino (prop_treino) / teste, ESTRATIFICADO por país.
#   - CV: PSU de treino repartidos em k_folds_cv folds, balanceados por país (para afinação).
#   - LOCO (leave-one-country-out): o país é o grupo de teste; treina nos outros 15.
#   - ESTRATO B (7 países ≥80%) entra inteiro como teste externo de transportabilidade.
#
# Chave: particoes.parquet liga-se a ds_clean por id_mulher. O desfecho e as covariáveis
# ficam em ds_clean; aqui guardam-se só as atribuições de partição (sem duplicar dados).
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(here)
})
source(here("2.code", "R", "config.R"))
source(here("2.code", "R", "funcoes.R"))

dc <- ler_parquet(file.path(path_process, "ds_clean.parquet")) |>
  filter(amostra == "principal") |>
  select(id_mulher, pais, nome_pais, estrato, psu_id, hiv_anc_resultado, w_norm)

verificar(sum(is.na(dc$hiv_anc_resultado)) == 0, "desfecho com NA na amostra principal")

set.seed(semente)

# ---------------------------------------------------------------------------
# 1) Holdout por PSU, estratificado por país (só Estrato A)
# ---------------------------------------------------------------------------
psu_A <- dc |> filter(estrato == "A") |> distinct(pais, psu_id)
psu_A <- psu_A |>
  group_by(pais) |>
  mutate(split = if_else(runif(n()) <= prop_treino, "treino", "teste")) |>
  ungroup()

# 2) Folds de CV (agrupados por PSU) dentro do treino, balanceados por país
psu_treino <- psu_A |> filter(split == "treino") |>
  group_by(pais) |>
  mutate(cv_fold = sample(rep_len(seq_len(k_folds_cv), n()))) |>
  ungroup() |>
  select(psu_id, cv_fold)

psu_A <- psu_A |> left_join(psu_treino, by = "psu_id")

# ---------------------------------------------------------------------------
# 3) Monta a tabela por mulher
# ---------------------------------------------------------------------------
part_A <- dc |> filter(estrato == "A") |>
  left_join(psu_A, by = c("pais", "psu_id")) |>
  mutate(particao = split)                       # "treino" / "teste"

part_B <- dc |> filter(estrato == "B") |>
  mutate(split = "externo_B", cv_fold = NA_integer_, particao = "externo_B")

particoes <- bind_rows(part_A, part_B) |>
  transmute(id_mulher, pais, nome_pais, estrato,
            psu_id, particao, cv_fold,
            loco_grupo = pais,                     # grupo de teste no LOCO
            desfecho = hiv_anc_resultado, w_norm)

# Sanidade: cada PSU num só lado; folds só no treino; nada em falta
chk_psu <- particoes |> filter(estrato == "A") |>
  group_by(psu_id) |> summarise(k = n_distinct(particao), .groups = "drop")
verificar(all(chk_psu$k == 1), "um PSU caiu em treino e teste (fuga por conglomerado)")
verificar(all(is.na(particoes$cv_fold[particoes$particao != "treino"])),
          "cv_fold atribuído fora do treino")
verificar(all(!is.na(particoes$cv_fold[particoes$particao == "treino"])),
          "mulher de treino sem cv_fold")

gravar_parquet(particoes, file.path(path_process, "particoes.parquet"))

# ---------------------------------------------------------------------------
# 4) Tabelas de resumo
# ---------------------------------------------------------------------------
resumo <- particoes |>
  group_by(nome_pais, estrato, particao) |>
  summarise(n = n(), psu = n_distinct(psu_id),
            prev = round(100 * weighted.mean(desfecho, w_norm), 1), .groups = "drop") |>
  arrange(estrato, nome_pais, particao)
write_csv(resumo, file.path(path_tabs, "05_particoes_resumo.csv"))

# Balanço global treino vs teste (Estrato A) — confirma prevalência preservada
balanco <- particoes |> filter(estrato == "A") |>
  group_by(particao) |>
  summarise(n = n(), psu = n_distinct(psu_id),
            prev = round(100 * weighted.mean(desfecho, w_norm), 1), .groups = "drop")

# LOCO: tamanho de cada fold (país = teste; treino = restantes A)
n_A <- sum(particoes$estrato == "A")
loco <- particoes |> filter(estrato == "A") |>
  group_by(loco_grupo, nome_pais) |>
  summarise(n_teste = n(), n_treino = n_A - n(),
            prev_teste = round(100 * weighted.mean(desfecho, w_norm), 1), .groups = "drop") |>
  arrange(nome_pais)
write_csv(loco, file.path(path_tabs, "05_loco_folds.csv"))

cat("Etapa 05 concluída.\n")
cat(sprintf("  Estrato A: %d mulheres, %d PSU, %d países.\n",
            n_A, n_distinct(particoes$psu_id[particoes$estrato == "A"]), 16L))
cat("  Balanço treino/teste (Estrato A):\n"); print(as.data.frame(balanco))
cat(sprintf("  Folds de CV: %d (agrupados por PSU, balanceados por país).\n", k_folds_cv))
cat(sprintf("  LOCO: 16 folds (1 país de teste, 15 de treino). Estrato B (%d mulheres) = teste externo.\n",
            sum(particoes$estrato == "B")))
cat("  Ficheiros: particoes.parquet, 05_particoes_resumo.csv, 05_loco_folds.csv\n")
