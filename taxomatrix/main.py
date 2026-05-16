from __future__ import annotations

import sys
import time
from pathlib import Path
import subprocess

import typer

from taxomatrix.export import export_label_mapping, export_matrix_csv, export_matrix_excel
from taxomatrix.fastani_runner import FastANIConfig, get_fastani_version, run_all_vs_all
from taxomatrix.parser import build_ani_matrix, parse_fastani_output
from taxomatrix.plotting import plot_clustered_heatmap
from taxomatrix.utils import (
    build_genome_label_mapping,
    discover_fasta_files,
    ensure_directory,
    format_command,
    validate_unique_display_labels,
    validate_unique_figure_labels,
    validate_unique_genome_labels,
    write_run_log,
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="TaxoMatrix: generate ANI matrices and clustered heatmaps from genome FASTA files.",
)


@app.callback()
def main() -> None:
    """TaxoMatrix command-line interface."""


@app.command()
def run(
    input_dir: Path = typer.Argument(..., exists=True, file_okay=False, dir_okay=True, resolve_path=True),
    output_dir: Path = typer.Option(
        Path("taxomatrix_output"),
        "--output-dir",
        "-o",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Directory for all TaxoMatrix outputs.",
    ),
    fastani_binary: str = typer.Option(
        "fastANI",
        "--fastani-binary",
        help="FastANI executable name or absolute path.",
    ),
    threads: int = typer.Option(1, "--threads", "-t", min=1, help="Threads to pass to FastANI."),
    frag_len: int | None = typer.Option(
        None,
        "--frag-len",
        min=1,
        help="Optional FastANI fragment length override.",
    ),
    min_fraction: float | None = typer.Option(
        None,
        "--min-fraction",
        min=0.0,
        max=1.0,
        help="Optional FastANI minimum shared genome fraction.",
    ),
    show_values: bool = typer.Option(
        False,
        "--show-values/--no-show-values",
        help="Show ANI values inside heatmap cells.",
    ),
    species_boundary: float = typer.Option(
        95.0,
        "--species-boundary",
        help="ANI reference value used when boundary annotation is enabled.",
    ),
    annotate_boundary: bool = typer.Option(
        False,
        "--annotate-boundary/--no-annotate-boundary",
        help="Add a species-boundary note and colorbar marker to the heatmap.",
    ),
    metadata: Path | None = typer.Option(
        None,
        "--metadata",
        exists=True,
        file_okay=True,
        dir_okay=False,
        resolve_path=True,
        help="Optional metadata TSV with genome, habitat, source, and pH columns.",
    ),
) -> None:
    """Run the complete TaxoMatrix ANI workflow."""
    start_time = time.perf_counter()
    warnings: list[str] = []
    try:
        ensure_directory(output_dir)
        genome_files = discover_fasta_files(input_dir)
        validate_unique_genome_labels(genome_files)
        if len(genome_files) < 2:
            raise ValueError(
                "Fewer than 2 genome FASTA files were found. "
                "TaxoMatrix requires at least 2 genomes for all-vs-all ANI."
        )
        label_mapping = build_genome_label_mapping(genome_files)
        validate_unique_display_labels(label_mapping)
        validate_unique_figure_labels(label_mapping)
        path_to_label = {
            path.resolve(): label
            for path, label in zip(genome_files, label_mapping["display_label"], strict=True)
        }
        typer.echo(f"Discovered {len(genome_files)} genome FASTA files in {input_dir}")

        config = FastANIConfig(
            executable=fastani_binary,
            threads=threads,
            frag_len=frag_len,
            min_fraction=min_fraction,
        )
        fastani_version = get_fastani_version(fastani_binary)

        fastani_result = run_all_vs_all(genome_files=genome_files, output_dir=output_dir, config=config)
        typer.echo(f"FastANI results written to {fastani_result.raw_output_path}")

        results = parse_fastani_output(fastani_result.raw_output_path)
        ani_matrix = build_ani_matrix(
            results,
            genome_files=genome_files,
            path_to_label=path_to_label,
        )
        missing_pair_count = int(ani_matrix.isna().sum().sum() // 2)
        if missing_pair_count:
            warnings.append(
                f"{missing_pair_count} genome pair(s) had no ANI result and were left blank in the matrix."
            )

        csv_path = export_matrix_csv(ani_matrix, output_dir / "ani_matrix.csv")
        mapping_path = export_label_mapping(label_mapping, output_dir / "genome_label_mapping.tsv")
        heatmap_artifacts = plot_clustered_heatmap(
            ani_matrix,
            output_dir,
            label_mapping=label_mapping,
            show_values=show_values,
            species_boundary=species_boundary,
            annotate_boundary=annotate_boundary,
            metadata_path=metadata,
        )
        warnings.extend(heatmap_artifacts.warnings)
        excel_path = export_matrix_excel(
            ani_matrix,
            output_dir / "ani_matrix.xlsx",
            metadata_sheet=heatmap_artifacts.metadata_sheet,
        )
        runtime_seconds = time.perf_counter() - start_time
        output_paths = [
            str(csv_path),
            str(excel_path),
            str(mapping_path),
            *(str(path) for path in heatmap_artifacts.output_paths),
        ]
        if heatmap_artifacts.metadata_legend_path is not None:
            output_paths.append(str(heatmap_artifacts.metadata_legend_path))
        output_paths.append(str(log_path := output_dir / "taxomatrix.log"))
        log_path = write_run_log(
            log_path,
            command=format_command([sys.executable, "-m", "taxomatrix.main", *sys.argv[1:]]),
            fastani_command=format_command(fastani_result.command),
            runtime_seconds=runtime_seconds,
            genome_count=len(genome_files),
            fastani_version=fastani_version,
            genome_labels=label_mapping["display_label"].tolist(),
            output_paths=output_paths,
            warnings=warnings,
        )

        typer.echo(f"ANI matrix CSV written to {csv_path}")
        typer.echo(f"ANI matrix Excel written to {excel_path}")
        typer.echo(f"Genome label mapping written to {mapping_path}")
        if heatmap_artifacts.metadata_legend_path is not None:
            typer.echo(f"Metadata legend written to {heatmap_artifacts.metadata_legend_path}")
        for heatmap_path in heatmap_artifacts.output_paths:
            typer.echo(f"Clustered heatmap written to {heatmap_path}")
        typer.echo(f"Run log written to {log_path}")
        for warning in warnings:
            typer.echo(f"Warning: {warning}", err=True)
    except subprocess.CalledProcessError as error:
        stderr = error.stderr.strip() if error.stderr else "No stderr captured."
        typer.echo(f"FastANI command failed:\n{stderr}", err=True)
        raise typer.Exit(code=1) from error
    except (FileNotFoundError, NotADirectoryError, ValueError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(code=1) from error


if __name__ == "__main__":
    app()
