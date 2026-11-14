# ============================================================================
# P06 - Testagem de HIV no pré-natal na África Subsaariana
# Etapa 02 - Harmonização das covariáveis, desenho amostral e pesos
# Autor: Manuel Mahoche
# Entrada: 1.data/process/ds_raw.parquet, 1.data/ref/regioes_v024.csv
# Saída:   1.data/process/ds_clean.parquet   (banco analítico, uma linha por mulher)
#          3.output/tabs/02_missing_por_pais.csv
#          3.output/tabs/02_dicionario_variaveis.csv
#
# Os códigos das categorias foram conferidos nos dicionários das fases 7 e 8
# (MZIR81FL.DO e Recode Manual DHS-7); são idênticos para todas as variáveis usadas.
# Ausência estrutural (dados do parceiro em mulheres fora de união) vira categoria própria.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(forcats); library(here)
})
source(here("2.code", "R", "config.R"))
source(here("2.code", "R", "funcoes.R"))

ds  <- ler_parquet(file.path(path_process, "ds_raw.parquet"))
reg <- read_csv(file.path(path_ref, "regioes_v024.csv"), show_col_types = FALSE)

sim_nao <- function(cond) factor(if_else(cond, "Sim", "Não"), levels = c("Não", "Sim"))

ocup_grupo <- function(x) {
  x <- as.numeric(x)
  case_when(
    x == 0                       ~ "Não trabalha",
    x %in% c(1, 2)               ~ "Profissional/administrativo",
    x %in% c(4, 5)               ~ "Agrícola",
    x %in% c(3, 6, 7, 8, 9, 96)  ~ "Comércio/serviços/manual",
    TRUE                         ~ NA_character_
  )
}
niv_ocup <- c("Não trabalha", "Profissional/administrativo", "Agrícola", "Comércio/serviços/manual")

dc <- ds |>
  left_join(reg, by = c("pais", "v024")) |>
  mutate(
    # ------------------------- identificação e desenho ----------------------
    id_mulher  = paste(pais, caseid, sep = "_"),
    psu_id     = paste(pais, v021, sep = "_"),
    estrato_id = paste(pais, v022, sep = "_"),
    regiao_id  = paste(pais, v024, sep = "_"),
    regiao_nome = coalesce(regiao_nome, as.character(v024)),
    ano_inquerito = as.numeric(v007),

    # ------------------------- sociodemográficas ----------------------------
    idade = as.numeric(v012),
    idade_grp = cut(idade, c(14, 24, 34, 49), labels = c("15-24", "25-34", "35-49")),
    uniao = factor(if_else(v502 == 1, "Em união", "Fora de união"),
                   levels = c("Fora de união", "Em união")),
    residencia = factor(case_when(v025 == 1 ~ "Urbana", v025 == 2 ~ "Rural"),
                        levels = c("Urbana", "Rural")),
    chefe_sexo = factor(case_when(v151 == 1 ~ "Masculino", v151 == 2 ~ "Feminino"),
                        levels = c("Masculino", "Feminino")),
    riqueza = factor(v190, levels = 1:5,
                     labels = c("Q1 (mais pobre)", "Q2", "Q3", "Q4", "Q5 (mais rico)")),
    escolaridade = factor(case_when(v106 == 0 ~ "Nenhuma", v106 == 1 ~ "Primária",
                                    v106 %in% 2:3 ~ "Secundária ou mais"),
                          levels = c("Nenhuma", "Primária", "Secundária ou mais")),
    ocupacao = factor(ocup_grupo(v717), levels = niv_ocup),
    educ_parceiro = factor(case_when(
      v502 != 1      ~ "Sem parceiro",
      v701 == 0      ~ "Nenhuma",
      v701 == 1      ~ "Primária",
      v701 %in% 2:3  ~ "Secundária ou mais",
      TRUE           ~ NA_character_),
      levels = c("Nenhuma", "Primária", "Secundária ou mais", "Sem parceiro")),
    ocup_parceiro = factor(if_else(v502 != 1, "Sem parceiro", ocup_grupo(v705)),
                           levels = c(niv_ocup, "Sem parceiro")),

    # ------------------------- informação ------------------------------------
    media_semanal = {
      m <- pmax(na_dhs(v157, 9), na_dhs(v158, 9), na_dhs(v159, 9), na.rm = TRUE)
      factor(case_when(m >= 2 ~ "Sim", m %in% 0:1 ~ "Não"), levels = c("Não", "Sim"))
    },
    internet_12m = factor(case_when(v171a == 1 ~ "Sim", v171a %in% c(0, 2, 3) ~ "Não"),
                          levels = c("Não", "Sim")),
    # Indicador DHS: sabe que o risco de transmissão vertical se reduz com medicamentos.
    # Quem não ouviu falar de HIV ou não foi perguntada pelo salto do questionário entra como "Não".
    conhece_ptv = sim_nao(v824 %in% 1),
    # CORREÇÃO (2026-10-29): "nunca ouviu falar de HIV" vinha de v750, que é "ouviu falar de uma ITS".
    # A pergunta sobre HIV/AIDS é v751 (ausente nos DHS-8 que não a fizeram: aí assume-se que todas
    # ouviram falar; v825 e v857a foram perguntadas a todas). v750 não ausente = módulo de HIV aplicado.
    estigma = factor(case_when(
      v751 %in% 0                    ~ "Não ouviu falar de HIV",
      v825 %in% 0 | v857a %in% 0     ~ "Atitude discriminatória",
      !is.na(v750)                   ~ "Sem atitude discriminatória"),
      levels = c("Sem atitude discriminatória", "Atitude discriminatória", "Não ouviu falar de HIV")),

    # ------------------------- reprodutivas ----------------------------------
    paridade = cut(as.numeric(v201), c(0, 1, 4, Inf), labels = c("1", "2-4", "5+")),
    idade_1parto = cut(as.numeric(v212), c(0, 17, 24, Inf), labels = c("<18", "18-24", "25+")),
    intencao_gravidez = factor(case_when(m10_1 == 1 ~ "Naquele momento", m10_1 == 2 ~ "Mais tarde",
                                         m10_1 == 3 ~ "Não queria mais"),
                               levels = c("Naquele momento", "Mais tarde", "Não queria mais")),
    desejo_filhos = factor(case_when(v602 == 1 ~ "Quer mais", v602 == 2 ~ "Indecisa",
                                     v602 %in% 3:5 ~ "Não quer/não pode"),
                           levels = c("Quer mais", "Indecisa", "Não quer/não pode")),

    # ------------------------- acesso ----------------------------------------
    distancia_problema = factor(case_when(v467d == 1 ~ "Problema grande",
                                          v467d %in% c(0, 2) ~ "Não é grande problema"),
                                levels = c("Não é grande problema", "Problema grande")),
    seguro_saude = factor(case_when(v481 == 1 ~ "Sim", v481 == 0 ~ "Não"), levels = c("Não", "Sim")),
    decisao_saude = factor(case_when(
      v502 != 1           ~ "Sem parceiro",
      v743a %in% 1:2      ~ "Participa",
      v743a %in% 4:6      ~ "Não participa"),
      levels = c("Participa", "Não participa", "Sem parceiro")),

    # ------------------------- pós-exposição (serviço) -----------------------
    anc_visitas = na_dhs(m14_1, c(98, 99)),
    anc_qualquer = factor(case_when(anc_visitas >= 1 ~ "Sim", anc_visitas == 0 ~ "Não"), levels = c("Não", "Sim")),
    anc_4mais = factor(case_when(anc_visitas >= 4 ~ "Sim", anc_visitas < 4 ~ "Não"), levels = c("Não", "Sim")),
    anc_inicio_trim = factor(case_when(
      anc_visitas == 0                    ~ "Sem pré-natal",
      na_dhs(m13_1, c(98, 99)) <= 3       ~ "1º trimestre",
      na_dhs(m13_1, c(98, 99)) %in% 4:10  ~ "Depois do 1º trimestre"),
      levels = c("1º trimestre", "Depois do 1º trimestre", "Sem pré-natal")),
    parto_instituicao = factor(case_when(m15_1 >= 20 & m15_1 < 96 ~ "Sim",
                                         m15_1 < 20 | m15_1 == 96 ~ "Não"),
                               levels = c("Não", "Sim")),
    testada_alguma_vez = sim_nao(v781 %in% 1),

    # ------------------------- pesos -----------------------------------------
    w_dhs = as.numeric(v005) / 1e6
  )

# Pesos: normalizado (média 1 por país) e de igual contribuição por país dentro de cada amostra
dc <- dc |>
  group_by(pais) |>
  mutate(w_norm = w_dhs / mean(w_dhs)) |>
  group_by(amostra) |>
  mutate(n_amostra = n(), k_paises = n_distinct(pais)) |>
  group_by(pais) |>
  mutate(w_igual = w_norm * (first(n_amostra) / first(k_paises)) / n()) |>
  ungroup() |>
  select(-n_amostra, -k_paises)

vars_saida <- c(
  "id_mulher", "pais", "nome_pais", "regiao_ssa", "ficheiro", "fase", "ano_inquerito", "amostra", "estrato",
  "psu_id", "estrato_id", "regiao_id", "regiao_nome", "v001", "v021", "v022", "v024",
  "w_dhs", "w_norm", "w_igual",
  "hiv_anc_resultado", "hiv_anc_testada", "hiv_teste_no_anc", "hiv_teste_no_parto",
  "idade", vars_preditores, "decisao_saude",
  "anc_visitas", vars_pos_exposicao, "testada_alguma_vez"
)
dc <- select(dc, all_of(vars_saida))

# Sanidade
verificar(all(abs(tapply(dc$w_norm, dc$pais, mean) - 1) < 1e-8), "w_norm não tem média 1 por país")
verificar(!any(duplicated(dc$id_mulher)), "id_mulher duplicado")

# Ausência por variável e país (%)
miss <- dc |>
  group_by(pais) |>
  summarise(across(c(all_of(vars_preditores), all_of(vars_pos_exposicao)),
                   ~ round(100 * mean(is.na(.x)), 1)), n = n(), .groups = "drop")

dicionario <- tibble(
  variavel = c(vars_preditores, "decisao_saude", vars_pos_exposicao),
  bloco = c(rep("Sociodemográfico", length(vars_bloco_sociodem)), rep("Informação", length(vars_bloco_info)),
            rep("Reprodutivo", length(vars_bloco_reprod)), rep("Acesso", length(vars_bloco_acesso)),
            "Acesso", rep("Pós-exposição (serviço)", length(vars_pos_exposicao))),
  categorias = sapply(variavel, function(v) if (is.factor(dc[[v]])) paste(levels(dc[[v]]), collapse = " | ") else "texto")
)

gravar_parquet(dc, file.path(path_process, "ds_clean.parquet"))
write_csv(miss, file.path(path_tabs, "02_missing_por_pais.csv"))
write_csv(dicionario, file.path(path_tabs, "02_dicionario_variaveis.csv"))

cat("Etapa 02 concluída:", nrow(dc), "mulheres,", ncol(dc), "variáveis\n")
