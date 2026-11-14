# Correr o pipeline P06 (testagem de HIV no pré-natal, DHS África Subsaariana)

## 0. Pré-requisitos
- R (>= 4.3) com os pacotes: arrow, dplyr, tidyr, purrr, readr, forcats, here, jsonlite.
  ```r
  install.packages(c("arrow", "dplyr", "tidyr", "purrr", "readr", "forcats", "here", "jsonlite"))
  ```
- Os IR do DHS em Parquet, um ficheiro por inquérito (`XXIRnnDT.parquet`).
  Exemplo de pasta: `CAMINHO/PARA/Parquet_DHS/IR` (definir em `DHS_PARQUET_IR`).

## 1. Apontar para os dados
O caminho dos IR está em `2.code/R/config.R` (`path_ir`). Noutra máquina, em vez de editar o código,
defina a variável de ambiente antes de abrir o R, ou no `.Renviron`:
```
DHS_PARQUET_IR=/caminho/para/Parquet_DHS/IR
```

## 2. Ordem de execução
Abrir o projeto (`P06_testagem-hiv-prenatal-ssa-ml.Rproj`) e correr:
```r
source("2.code/00_executar_tudo.R")
```
ou etapa a etapa:
```
Rscript 2.code/00_acesso.R                   # ~1 min; grava 1.data/raw/dhs_ir_anc.parquet
Rscript 2.code/01_elegibilidade_desfecho.R   # pára se o desfecho divergir > 1 pp do STATcompiler
Rscript 2.code/02_limpeza.R                  # grava 1.data/process/ds_clean.parquet
```

## 3. Verificação (o que deve sair)
- `01_qc_desfecho_statcompiler.csv`: os 28 inquéritos com `dif_pp = 0` e `qc_ok = TRUE`.
- `ds_clean.parquet`: 116 961 mulheres (99 285 na amostra principal, 23 países;
  17 676 na amostra de sensibilidade, 5 países), 51 variáveis.
- `mean(w_norm) = 1` em cada país; soma de `w_igual` igual entre países da mesma amostra.

## 4. Acrescentar ou trocar um inquérito
1. Pôr o novo `XXIRnnDT.parquet` na pasta dos IR.
2. Editar a linha do país em `inqueritos` (`config.R`) e o valor do STATcompiler em `statcompiler`.
3. Acrescentar os nomes das regiões em `1.data/ref/regioes_v024.csv` (rótulos de `v024` do `.DTA`).
4. Correr de novo as etapas 00 a 02. O QC confirma se o desfecho bate com o valor oficial.
