# TaxoMatrix

TaxoMatrix is a Python command-line workflow for generating all-vs-all ANI outputs from a folder of microbial genome FASTA files. It combines FastANI execution, publication-friendly genome labeling, matrix export, and clustered heatmap rendering in a reproducible local workflow.

## Features

- all-vs-all ANI using FastANI
- publication-friendly genome labels with a separate short figure label
- ANI matrix export to CSV and Excel
- clustered heatmap export to PNG, SVG, and PDF
- optional metadata sidebars for `habitat`, `source`, and `pH`
- run logging and output provenance
- modular code structure that can be extended later for AAI

## Project structure

```text
taxomatrix/
  __init__.py
  main.py
  fastani_runner.py
  parser.py
  plotting.py
  export.py
  utils.py
LICENSE
CITATION.cff
environment.yml
example_metadata.tsv
EXAMPLE_COMMAND_OUTPUTS.md
requirements.txt
README.md
```

## Requirements

- Python 3.11
- FastANI installed and available on `PATH`, or provided explicitly with `--fastani-binary`

Create an environment with Conda:

```bash
conda env create -f environment.yml
conda activate taxomatrix
```

Install Python dependencies:

```bash
python3.11 -m pip install -r requirements.txt
```

## Usage

Basic run:

```bash
python3.11 -m taxomatrix.main run /path/to/genome_fastas
```

Write outputs to a custom directory:

```bash
python3.11 -m taxomatrix.main run /path/to/genome_fastas --output-dir results/taxomatrix_run
```

Use a specific FastANI binary and more threads:

```bash
python3.11 -m taxomatrix.main run /path/to/genome_fastas \
  --fastani-binary /usr/local/bin/fastANI \
  --threads 8
```

Optional FastANI parameters:

```bash
python3.11 -m taxomatrix.main run /path/to/genome_fastas \
  --frag-len 3000 \
  --min-fraction 0.2
```

Show ANI values in heatmap cells:

```bash
python3.11 -m taxomatrix.main run /path/to/genome_fastas \
  --output-dir results/taxomatrix_values \
  --show-values
```

Annotate the species boundary note:

```bash
python3.11 -m taxomatrix.main run /path/to/genome_fastas \
  --annotate-boundary \
  --species-boundary 95.0
```

Use metadata sidebars:

```bash
python3.11 -m taxomatrix.main run /path/to/genome_fastas \
  --metadata example_metadata.tsv
```

## Outputs

Each run writes:

- `fastani_raw.tsv`
- `ani_matrix.csv`
- `ani_matrix.xlsx`
- `ani_clustered_heatmap.png`
- `ani_clustered_heatmap.svg`
- `ani_clustered_heatmap.pdf`
- `genome_label_mapping.tsv`
- `taxomatrix.log`
- `metadata_legend.tsv` when `--metadata` is supplied

## Example metadata

An example metadata file is included as [example_metadata.tsv](/home/justus/example_metadata.tsv:1):

```tsv
genome	habitat	source	pH
GCA_000188155.3	soil	culture	6.8
GCA_000283235.1	wetland	isolate	6.3
GCA_000304315.1	peatland	enrichment	5.9
GCA_000372845.1	freshwater	isolate	7.1
```

## Example command outputs

Representative command outputs are documented in [EXAMPLE_COMMAND_OUTPUTS.md](/home/justus/EXAMPLE_COMMAND_OUTPUTS.md:1).

## Notes

- Missing ANI values are left blank in the exported matrix and treated as zero similarity for clustering.
- Output tables use publication-friendly genome labels derived from FASTA metadata and assembly accessions.
- Figure axes use abbreviated scientific labels while full labels remain in tables and logs.
- The Excel workbook includes an `ANI Matrix` sheet and a `Summary` sheet with off-diagonal ANI statistics, plus a `Metadata` sheet when metadata is supplied.
- The code is organized so an AAI runner and parser can be added alongside the ANI workflow later.

## Version

This release is intended as `TaxoMatrix v0.1.0`.
