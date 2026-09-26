# PPI_enrichment_interactif

Interactive and command-line workflow for high-confidence STRING PPI contextualization and g:Profiler functional over-representation analysis.

## Paper implementation

The repository packages `PPI_enrichment_single_final_v1_0.py` (v1.0.0), the PPI/enrichment script supplied for the manuscript **Evidence-guided target prioritization in pulmonary thrombo-inflammatory disease**. In the reported study, the fixed shortlist was `REN, SELP, ANGPT2, CXCL10, IFNG`; primary ORA used the fixed study-derived 197-target background; STRING used Homo sapiens taxonomy 9606, required score 700, and up to 5 first-shell partners per seed. Network expansion is descriptive context and does not re-rank targets.

**Important:** the repository does not contain the study's real 197-target matrix unless you add/deposit it. To reproduce the publication, use the exact frozen `target_evidence_matrix.csv` from the target-discovery output/reproducibility package. The included example CSV is only a format example and must not be used to reproduce manuscript results.

## What is interactive?

Run the script without the corresponding command-line values and it prompts for: (1) shortlisted targets, (2) target-discovery output root, and (3) output path. The target set is therefore replaceable. The background can be supplied with `--background`; STRING score and first-shell size are configurable CLI parameters.

The implementation is reusable for other **human** gene/protein sets. It is not presently species-generic because STRING taxonomy `9606` and g:Profiler organism `hsapiens` are constants in v1.0.0. Reuse in another disease is possible with a suitable disease/study-specific background. Cross-species reuse requires code modification and validation.

## Requirements

- Python 3.9+ recommended
- Internet access during execution (STRING and g:Profiler APIs)
- Python packages in `requirements.txt`
- A background CSV with one of these columns: `Target`, `HGNC_Symbol`, `Gene`, or `Symbol`

## Install

### Linux/macOS
```bash
git clone https://github.com/Salma2306/PPI_enrichment_interactif.git
cd PPI_enrichment_interactif
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

### Windows PowerShell
```powershell
git clone https://github.com/Salma2306/PPI_enrichment_interactif.git
cd PPI_enrichment_interactif
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Reproduce the paper PPI/enrichment step

```bash
python PPI_enrichment_single_final_v1_0.py \
  --targets "REN;SELP;ANGPT2;CXCL10;IFNG" \
  --background /absolute/path/to/target_evidence_matrix.csv \
  --out outputs/paper_run \
  --string-score 700 \
  --partners-per-seed 5
```

For the publication run, the supplied background should contain exactly 197 study-derived targets. The script warns if another background size is supplied, but intentionally does not alter the input.

## Interactive run

```bash
python PPI_enrichment_single_final_v1_0.py
```

Example prompts:
```text
Shortlisted targets [REN; SELP; ANGPT2; CXCL10; IFNG]:
Target-discovery output root [...]:
```

If `--background` is omitted, the script searches the output root for `06_cross_database/target_evidence_matrix.csv` or `target_evidence_matrix.csv`.

## Run on another human target set

```bash
python PPI_enrichment_single_final_v1_0.py \
  --targets "IL6;TNF;VWF;F2" \
  --background /absolute/path/to/your_study_background.csv \
  --out outputs/example_custom \
  --string-score 700 \
  --partners-per-seed 5
```

## Outputs

The workflow produces:

- `01_STRING_mapping.csv`
- `02_seed_PPI_edges.csv`
- `03_seed_enrichment_vs_197.csv` (filename retained from the paper implementation)
- `04_seed_enrichment_FDR_significant.csv`
- `05_first_shell_context.csv`
- `06_expanded_PPI_edges.csv`
- `07_expanded_enrichment_context_only.csv`
- `08_background_universe.csv`
- `09_QC_provenance.json`
- `10_gProfiler_requests.json`
- `Figure_3_PPI_enrichment.png`
- `Figure_3_PPI_enrichment.pdf`
- `Figure_3_PPI_enrichment.svg`

`09_QC_provenance.json` records the software version, requested/mapped seeds, background file and size, STRING threshold, first-shell setting, enrichment count and network dimensions. `10_gProfiler_requests.json` preserves the API request bodies.

## Statistical safeguards

- Primary ORA uses the supplied study-derived custom background.
- Multiple testing uses FDR with threshold 0.05.
- All corrected results are retained.
- A negative primary ORA is retained rather than forcing significance.
- First-shell STRING expansion is contextual only.
- No permutation, leave-one-out target reranking, adaptive thresholding or forced significance is performed by this script.

## External services

- STRING API: version-specific v12.0 endpoint is used by v1.0.0.
- g:Profiler g:GOSt API: `hsapiens`, custom annotated background for primary ORA, FDR correction.

Because g:Profiler annotations are an external live resource, exact long-term replication should archive the returned CSV/JSON outputs and execution date with the publication release.

## Quick validation

```bash
python -m py_compile PPI_enrichment_single_final_v1_0.py
python PPI_enrichment_single_final_v1_0.py --help
```

## Citation

Manuscript citation and Zenodo DOI should be added once available. For a publication release, create a tagged GitHub release (for example `v1.0.0`) and archive that release with Zenodo.

## License

MIT.
