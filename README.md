# Antenatal HIV testing in sub-Saharan Africa: analysis code (DHS)

Code for a multi-country analysis of HIV testing with receipt of results during antenatal care or delivery,
using Demographic and Health Surveys (DHS) Individual Recode files from sub-Saharan Africa.
The analysis develops machine learning prediction models, evaluates how well they transport across
countries (leave-one-country-out and external validation), and measures socioeconomic inequality.

> **Status.** This repository contains the code for the main prediction and inequality pipeline
> (stages 00-07, 10, 11, 15, 18, 19, 24). Code for additional analyses will be added when the
> corresponding manuscripts are published. Results are reported in the associated article
> (citation and DOI to be added after publication).

## Data are not included

DHS microdata cannot be redistributed. To reproduce the analysis you need to:

1. Register at <https://dhsprogram.com> and request access to the Individual Recode (IR) files
   of the surveys listed in `2.code/R/config.R` (table `inqueritos`).
2. Convert each IR file to Parquet (one file per survey, named `XXIRnnDT.parquet`).
3. Point the code to that folder with the environment variable `DHS_PARQUET_IR`
   (or edit `path_ir` in `2.code/R/config.R`).

## Requirements

- R (4.x) and the packages `arrow`, `dplyr`, `tidyr`, `purrr`, `readr`, `forcats`, `tibble`, `here`,
  `jsonlite`, `survey`, `srvyr`, `glmnet`, `ranger`, `xgboost`, `pROC`, `ggplot2`, `scales`.
- Enough disk space for the stacked Parquet file written to `1.data/raw/` (not versioned).

## How to run

See `COMO_RODAR_OUTRA_MAQUINA.md` (in Portuguese) for step-by-step instructions. In short:

```bash
export DHS_PARQUET_IR=/path/to/Parquet_DHS/IR
bash 2.code/00_rodar_tudo.sh          # runs the stages in numeric order, one log per stage
```

Outputs are written to `3.output/` (tables, figures, models), which is not versioned.

## Pipeline stages

| Stage | Script | Purpose |
|---|---|---|
| 00 | `00_acesso.R` | Read the IR files and stack the study variables |
| 01 | `01_elegibilidade_desfecho.R` | Eligible population, outcome, and check against the STATcompiler indicator |
| 02 | `02_limpeza.R` | Harmonise covariates; survey design and weights |
| 03 | `03_descritiva.R` | Design-based descriptive statistics |
| 04 | `04_desigualdades.R` | Inequality measures by wealth |
| 05 | `05_particoes.R` | Validation splits (by primary sampling unit, cross-validation, leave-one-country-out) |
| 06 | `06_modelos.R` | Models and transportability |
| 07 | `07_calibracao_interpretabilidade.R` | Calibration and interpretability |
| 10 | `10_correlacao_preditores.R` | Predictor correlation and collinearity |
| 11 | `11_sensibilidade_dados.R` | Sensitivity to missing data and weighting |
| 15 | `15_sensibilidade_recalibracao.R` | Recalibration sensitivity |
| 18 | `18_desigualdade_corrigida.R` | Corrected concentration indices (Wagstaff, Erreygers) |
| 19 | `19_revisao_modelos.R` | Model uncertainty, equity and clinical utility |
| 24 | `24_sens_nunca_ouviu.R` | Analysis restricted to women who had heard of HIV |

Configuration (countries, predictors, seeds, stratum threshold) is in `2.code/R/config.R`.
`CORTE_ESTRATO` (environment variable) changes the coverage threshold that defines the strata.

## Notes

- Script comments and the run guide are in Portuguese.
- This is research code, not a software package; it is provided as is.
- Licence: to be defined (suggestion: MIT for code).
- Contact: to be added.
