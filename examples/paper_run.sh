#!/usr/bin/env bash
set -e
python PPI_enrichment_single_final_v1_0.py \
  --targets "REN;SELP;ANGPT2;CXCL10;IFNG" \
  --background /ABSOLUTE/PATH/TO/target_evidence_matrix.csv \
  --out outputs/paper_run \
  --string-score 700 \
  --partners-per-seed 5
