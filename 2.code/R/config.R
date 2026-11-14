# ============================================================================
# P06 - Testagem de HIV no pré-natal na África Subsaariana (ML + cenários)
# Configuração única do pipeline: caminhos, inquéritos, variáveis e referências
# Uso: source(here::here("2.code", "R", "config.R")) no início de cada etapa
# ============================================================================

# ----------------------------- Caminhos --------------------------------------
# Pasta com os IR em Parquet (um ficheiro por inquérito: XXIRnnDT.parquet).
# Pode ser redefinida sem editar o código com a variável de ambiente DHS_PARQUET_IR.
path_ir <- Sys.getenv(
  "DHS_PARQUET_IR",
  unset = "CAMINHO/PARA/Parquet_DHS/IR"   # ajustar, ou definir a variável de ambiente DHS_PARQUET_IR
)

path_raw     <- here::here("1.data", "raw")
path_process <- here::here("1.data", "process")
path_ref     <- here::here("1.data", "ref")
path_tabs    <- here::here("3.output", "tabs")
path_figs    <- here::here("3.output", "figs")

# ----------------------------- Inquéritos ------------------------------------
# DHS padrão mais recente por país com módulo de testagem de HIV no pré-natal.
# amostra: "principal" = trabalho de campo iniciado em 2018 ou depois;
#          "sensibilidade" = janela ampliada (2015-2017), só para análises de sensibilidade.
# Excluídos (documentado no README): Etiópia (ETIR81 = Mini-DHS 2019, sem módulo de HIV),
# Chade (TDIR71, 2014-15, módulo quase vazio), e países sem DHS desde 2015.
inqueritos <- tibble::tribble(
  ~pais, ~nome_pais,          ~ficheiro,   ~anos,       ~regiao_ssa, ~amostra,
  "AO",  "Angola",            "AOIR81DT",  "2023-24",   "Central",   "principal",
  "BF",  "Burkina Faso",      "BFIR81DT",  "2021",      "Ocidental", "principal",
  "CD",  "RD Congo",          "CDIR81DT",  "2023-24",   "Central",   "principal",
  "CI",  "Costa do Marfim",   "CIIR81DT",  "2021",      "Ocidental", "principal",
  "CM",  "Camarões",          "CMIR71DT",  "2018",      "Central",   "principal",
  "GA",  "Gabão",             "GAIR71DT",  "2019-21",   "Central",   "principal",
  "GH",  "Gana",              "GHIR8CDT",  "2022",      "Ocidental", "principal",
  "GM",  "Gâmbia",            "GMIR81DT",  "2019-20",   "Ocidental", "principal",
  "GN",  "Guiné",             "GNIR71DT",  "2018",      "Ocidental", "principal",
  "KE",  "Quénia",            "KEIR8CDT",  "2022",      "Oriental",  "principal",
  "LB",  "Libéria",           "LBIR7ADT",  "2019-20",   "Ocidental", "principal",
  "LS",  "Lesoto",            "LSIR81DT",  "2023-24",   "Austral",   "principal",
  "MD",  "Madagáscar",        "MDIR81DT",  "2021",      "Oriental",  "principal",
  "ML",  "Mali",              "MLIR8ADT",  "2023-24",   "Ocidental", "principal",
  "MR",  "Mauritânia",        "MRIR71DT",  "2019-21",   "Ocidental", "principal",
  "MW",  "Malawi",            "MWIR81DT",  "2024",      "Austral",   "principal",
  "MZ",  "Moçambique",        "MZIR81DT",  "2022-23",   "Austral",   "principal",
  "NG",  "Nigéria",           "NGIR8BDT",  "2023-24",   "Ocidental", "principal",
  "RW",  "Ruanda",            "RWIR81DT",  "2019-20",   "Oriental",  "principal",
  "SL",  "Serra Leoa",        "SLIR7ADT",  "2019",      "Ocidental", "principal",
  "SN",  "Senegal",           "SNIR8SDT",  "2023",      "Ocidental", "principal",
  "TZ",  "Tanzânia",          "TZIR82DT",  "2022",      "Oriental",  "principal",
  "ZM",  "Zâmbia",            "ZMIR81DT",  "2024",      "Austral",   "principal",
  "BJ",  "Benim",             "BJIR71DT",  "2017-18",   "Ocidental", "sensibilidade",
  "BU",  "Burundi",           "BUIR71DT",  "2016-17",   "Oriental",  "sensibilidade",
  "UG",  "Uganda",            "UGIR7BDT",  "2016",      "Oriental",  "sensibilidade",
  "ZA",  "África do Sul",     "ZAIR71DT",  "2016",      "Austral",   "sensibilidade",
  "ZW",  "Zimbabué",          "ZWIR72DT",  "2015",      "Austral",   "sensibilidade"
)

# Valores oficiais do STATcompiler (indicador HA_CATH_W_ATR: testadas para HIV no pré-natal
# ou no parto e que receberam o resultado, mulheres com nascido vivo nos 2 anos anteriores).
# Consultados em 16/10/2026. Servem de controlo de qualidade do desfecho (Etapa 01).
statcompiler <- tibble::tribble(
  ~pais, ~stat_pct, ~stat_n,
  "AO", 49.9,  4551, "BF", 54.9,  4721, "CD", 26.4,  9307, "CI", 50.8,  4184,
  "CM", 75.3,  3788, "GA", 89.0,  1655, "GH", 72.4,  3850, "GM", 58.1,  3441,
  "GN", 22.7,  3066, "KE", 94.8,  3985, "LB", 63.1,  2284, "LS", 91.3,  1087,
  "MD", 16.4,  4946, "ML", 26.6,  5916, "MR", 10.2,  4465, "MW", 96.5,  4841,
  "MZ", 62.7,  3652, "NG", 33.7, 10522, "RW", 96.7,  3135, "SL", 61.4,  3970,
  "SN", 43.3,  4071, "TZ", 88.2,  4308, "ZM", 90.7,  3540, "BJ", 40.2,  2680,
  "BU", 89.7,  5261, "UG", 91.5,  5992, "ZA", 94.5,  1376, "ZW", 90.3,  2367
)
tolerancia_qc_pp <- 1.0   # diferença máxima aceite (pontos percentuais)

# Corte dos estratos (decisão de desenho, fixada antes da modelagem)
corte_estrato_pct <- as.numeric(Sys.getenv("CORTE_ESTRATO", "80"))   # A: cobertura < corte (principal: 80); sensibilidade do Plano: 70 e 90

# Partições de ML (Etapa 05): semente única e parâmetros do desenho de validação.
semente       <- 20261020  # reprodutibilidade de todas as partições e modelos
prop_treino   <- 0.75      # fração de PSU para treino no holdout bloqueado por PSU
k_folds_cv    <- 5         # folds de validação cruzada (agrupados por PSU) dentro do treino

# ----------------------------- Variáveis a ler do IR -------------------------
vars_desenho <- c("caseid", "v000", "v001", "v002", "v003", "v005", "v007", "v008",
                  "v021", "v022", "v023", "v024", "v025")
# v750 = "ever heard of an STI" (usada só como indicador de que o módulo de HIV foi aplicado);
# v751 = "ever heard of AIDS" (ausente nos DHS-8 sem a pergunta: assume-se que todas ouviram falar).
vars_elegib  <- c("b3_01", "b19_01", "v208", "v750", "v751", "scamps", "sshort")
vars_desfecho <- c("v840", "v840a", "v841", "v841a", "m14_1")
vars_covar <- c(
  "v012", "v106", "v190", "v151", "v501", "v502", "v701", "v705", "v717",  # sociodemográficas
  "v157", "v158", "v159", "v171a",                                         # meios / internet
  "v824", "v825", "v857a",                                                 # conhecimento PTV / estigma
  "v201", "v212", "v602", "m10_1", "m13_1", "m15_1",                       # reprodutivas / serviço
  "v467d", "v481", "v743a", "v781"                                         # acesso / autonomia
)
vars_ler <- unique(c(vars_desenho, vars_elegib, vars_desfecho, vars_covar))

# Variáveis de texto (string). O formato varia entre países (umas vêm como
# character, outras como numérico/rotulado), pelo que são forçadas a character
# antes de empilhar os inquéritos, para evitar conflitos de classe no rbind.
vars_texto <- c("caseid", "v000", "scamps", "sshort")

# ----------------------------- Dicionário analítico --------------------------
# Preditores candidatos (os 16 do estudo de Burkina Faso, adaptados aos itens comuns às fases 7 e 8)
vars_bloco_sociodem <- c("idade_grp", "uniao", "residencia", "regiao_id", "chefe_sexo",
                         "riqueza", "escolaridade", "ocupacao", "educ_parceiro", "ocup_parceiro")
vars_bloco_info     <- c("media_semanal", "internet_12m", "conhece_ptv", "estigma")
vars_bloco_reprod   <- c("paridade", "idade_1parto", "intencao_gravidez", "desejo_filhos")
vars_bloco_acesso   <- c("distancia_problema", "seguro_saude")
# Pós-exposição no DAG (descritivas / mediadoras; não entram como preditores do modelo causal)
vars_pos_exposicao  <- c("anc_qualquer", "anc_4mais", "anc_inicio_trim", "parto_instituicao")

vars_preditores <- c(vars_bloco_sociodem, vars_bloco_info, vars_bloco_reprod, vars_bloco_acesso)
