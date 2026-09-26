#!/usr/bin/env bash
set -e
python PPI_enrichment_single_final_v1_0.py \
  --targets "GENE1;GENE2;GENE3" \
  --background /ABSOLUTE/PATH/TO/your_background.csv \
  --out outputs/custom_run \
  --string-score 700 \
  --partners-per-seed 5
