#!/bin/bash
# Corre as etapas 00-20 por ordem, com log por etapa e paragem na primeira falha.
# Uso: DHS_PARQUET_IR=/caminho/dos/IR bash 2.code/00_rodar_tudo.sh [etapa_inicial]
cd "$(dirname "$0")/.." || exit 1
inicio=${1:-00}
mkdir -p 3.output/logs
for f in 2.code/[0-2][0-9]_*.R; do
  n=$(basename "$f" | cut -c1-2)
  [[ "$n" < "$inicio" ]] && continue
  [[ "$(basename "$f")" == "00_executar_tudo.R" || "$f" == *_en.R ]] && continue
  t0=$(date +%s)
  # uma repetição: no OneDrive, falhas transitórias de leitura/gravação acontecem
  if Rscript "$f" > "3.output/logs/$(basename "$f" .R).log" 2>&1 || Rscript "$f" > "3.output/logs/$(basename "$f" .R).log" 2>&1; then
    echo "OK   $(basename "$f") ($(( $(date +%s) - t0 ))s)"
  else
    echo "FAIL $(basename "$f") -- ver 3.output/logs/$(basename "$f" .R).log"; exit 1
  fi
done
echo "TUDO CONCLUIDO"
