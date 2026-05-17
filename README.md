# TaxoMatrix

[![GitHub release](https://img.shields.io/github/v/release/justusnweze/TaxoMatrix)](https://github.com/justusnweze/TaxoMatrix/releases)

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
  downloader.py
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
- NCBI datasets CLI for accession-download mode (`--accessions`)

Create an environment with Conda:

```bash
conda env create -f environment.yml
conda activate taxomatrix
```

Install Python dependencies:

```bash
python3.11 -m pip install -r requirements.txt
```

Install the NCBI datasets CLI if you want to download genomes by accession:

```bash
conda install -c conda-forge ncbi-datasets-cli
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

Run from an accession file:

```bash
python3.11 -m taxomatrix.main run \
  --accessions genomes.txt \
  --output-dir results/taxomatrix_accessions
```

Combine local genomes with downloaded accessions:

```bash
python3.11 -m taxomatrix.main run /path/to/genome_fastas \
  --accessions genomes.txt \
  --output-dir results/taxomatrix_combined
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
- `downloaded_genomes/` when `--accessions` is supplied
- `CITATION_AND_ACKNOWLEDGEMENT.txt`

## Example metadata

An example metadata file is included as [example_metadata.tsv](example_metadata.tsv):

```tsv
genome	habitat	source	pH
GCA_000188155.3	soil	culture	6.8
GCA_000283235.1	wetland	isolate	6.3
GCA_000304315.1	peatland	enrichment	5.9
GCA_000372845.1	freshwater	isolate	7.1
```

## Example command outputs

Representative command outputs are documented in [EXAMPLE_COMMAND_OUTPUTS.md](EXAMPLE_COMMAND_OUTPUTS.md).

## Accession mode

TaxoMatrix can also download public assemblies from NCBI using an accession file containing `GCA_` or `GCF_` assembly accessions, one per line. Downloaded genomes are cached in `output_dir/downloaded_genomes/` so reruns do not redownload successful assemblies unnecessarily.

## Dependencies and Citation

TaxoMatrix builds upon established open-source bioinformatics and scientific Python tools. For transparent and reproducible scientific reporting, please cite both TaxoMatrix and the upstream software used in your workflow.

### Please also cite upstream tools

- FastANI
  Used for pairwise Average Nucleotide Identity calculations.
- pandas
  Used for matrix handling and tabular data processing.
- matplotlib
  Used for figure rendering.
- seaborn
  Used for clustered heatmap visualization.
- scipy
  Used for hierarchical clustering.
- openpyxl
  Used for Excel workbook export.
- NCBI datasets CLI (if accession mode is used)
  Used for downloading public genome assemblies.

TaxoMatrix does not replace citation of the underlying computational tools used in the workflow.

Detailed dependency citation placeholders are collected in [DEPENDENCY_CITATIONS.md](DEPENDENCY_CITATIONS.md).

## Citation and DOI

TaxoMatrix can be cited using the metadata in [CITATION.cff](CITATION.cff). A Zenodo DOI can be created later from GitHub releases to make the release record citable and archived.

## Notes

- Missing ANI values are left blank in the exported matrix and treated as zero similarity for clustering.
- Output tables use publication-friendly genome labels derived from FASTA metadata and assembly accessions.
- Figure axes use abbreviated scientific labels while full labels remain in tables and logs.
- The Excel workbook includes an `ANI Matrix` sheet and a `Summary` sheet with off-diagonal ANI statistics, plus a `Metadata` sheet when metadata is supplied.
- The code is organized so an AAI runner and parser can be added alongside the ANI workflow later.
- ANI results should be interpreted alongside genome quality, phylogenomics, phenotype, and nomenclatural rules.

## Version

This release is intended as `TaxoMatrix v0.1.1`.

### v0.1.1

TaxoMatrix v0.1.1 adds accession-based genome downloads, upstream citation guidance, and release metadata cleanup while preserving the original ANI workflow.
