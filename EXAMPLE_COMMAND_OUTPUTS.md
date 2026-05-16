# Example Command Outputs

## Baseline run

Command:

```bash
python -m taxomatrix.main run test_data/genomes -o test_output_refined -t 14
```

Example console output:

```text
Discovered 4 genome FASTA files in test_data/genomes
FastANI results written to test_output_refined/fastani_raw.tsv
ANI matrix CSV written to test_output_refined/ani_matrix.csv
ANI matrix Excel written to test_output_refined/ani_matrix.xlsx
Genome label mapping written to test_output_refined/genome_label_mapping.tsv
Clustered heatmap written to test_output_refined/ani_clustered_heatmap.png
Clustered heatmap written to test_output_refined/ani_clustered_heatmap.svg
Clustered heatmap written to test_output_refined/ani_clustered_heatmap.pdf
Run log written to test_output_refined/taxomatrix.log
```

## Heatmap values and species boundary

Command:

```bash
python -m taxomatrix.main run test_data/genomes -o test_output_values -t 14 --show-values --annotate-boundary
```

Example console output:

```text
Discovered 4 genome FASTA files in test_data/genomes
FastANI results written to test_output_values/fastani_raw.tsv
ANI matrix CSV written to test_output_values/ani_matrix.csv
ANI matrix Excel written to test_output_values/ani_matrix.xlsx
Genome label mapping written to test_output_values/genome_label_mapping.tsv
Clustered heatmap written to test_output_values/ani_clustered_heatmap.png
Clustered heatmap written to test_output_values/ani_clustered_heatmap.svg
Clustered heatmap written to test_output_values/ani_clustered_heatmap.pdf
Run log written to test_output_values/taxomatrix.log
```

## Metadata-enabled run

Command:

```bash
python -m taxomatrix.main run test_data/genomes -o test_output_metadata -t 14 --metadata test_data/metadata.tsv
```

Example console output:

```text
Discovered 4 genome FASTA files in test_data/genomes
FastANI results written to test_output_metadata/fastani_raw.tsv
ANI matrix CSV written to test_output_metadata/ani_matrix.csv
ANI matrix Excel written to test_output_metadata/ani_matrix.xlsx
Genome label mapping written to test_output_metadata/genome_label_mapping.tsv
Metadata legend written to test_output_metadata/metadata_legend.tsv
Clustered heatmap written to test_output_metadata/ani_clustered_heatmap.png
Clustered heatmap written to test_output_metadata/ani_clustered_heatmap.svg
Clustered heatmap written to test_output_metadata/ani_clustered_heatmap.pdf
Run log written to test_output_metadata/taxomatrix.log
```
