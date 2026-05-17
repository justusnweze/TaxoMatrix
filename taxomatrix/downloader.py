from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import zipfile

from taxomatrix.utils import ensure_directory, run_command


@dataclass(slots=True)
class DownloadArtifacts:
    genome_files: list[Path]
    warnings: list[str]
    notes: list[str]


def _find_genome_fasta(extract_dir: Path) -> Path:
    fasta_files = sorted(extract_dir.rglob("*_genomic.fna"))
    if not fasta_files:
        fasta_files = sorted(extract_dir.rglob("*.fna"))
    if not fasta_files:
        raise FileNotFoundError(
            f"Downloaded dataset did not contain a genomic FASTA file in {extract_dir}."
        )
    return fasta_files[0]


def download_accession_genomes(
    accessions: list[str],
    output_dir: Path,
    *,
    datasets_binary: str = "datasets",
) -> DownloadArtifacts:
    datasets_candidate = shutil.which(datasets_binary)
    if datasets_candidate:
        datasets_executable = datasets_candidate
    else:
        candidate_path = Path(datasets_binary)
        if candidate_path.exists() and candidate_path.is_file():
            datasets_executable = str(candidate_path)
        else:
            raise FileNotFoundError(
                f"NCBI datasets CLI executable not found: {datasets_binary}."
            )
    download_dir = ensure_directory(output_dir / "downloaded_genomes")
    warnings: list[str] = []
    notes: list[str] = []
    genome_files: list[Path] = []

    for accession in accessions:
        final_fasta_path = download_dir / f"{accession}_genomic.fna"
        if final_fasta_path.exists():
            genome_files.append(final_fasta_path)
            notes.append(f"Used cached downloaded genome for accession {accession}: {final_fasta_path}")
            continue

        accession_dir = ensure_directory(download_dir / accession)
        zip_path = accession_dir / f"{accession}.zip"
        extract_dir = accession_dir / "extract"
        if extract_dir.exists():
            shutil.rmtree(extract_dir)
        ensure_directory(extract_dir)

        command = [
            datasets_executable,
            "download",
            "genome",
            "accession",
            accession,
            "--include",
            "genome",
            "--filename",
            str(zip_path),
            "--no-progressbar",
        ]

        try:
            run_command(command, cwd=output_dir)
            with zipfile.ZipFile(zip_path) as archive:
                archive.extractall(extract_dir)
            genome_fasta = _find_genome_fasta(extract_dir)
            shutil.copy2(genome_fasta, final_fasta_path)
            genome_files.append(final_fasta_path)
            notes.append(f"Downloaded accession {accession} to {final_fasta_path}")
        except Exception as error:
            warnings.append(
                f"Accession download failed for {accession}: {error}"
            )
        finally:
            if zip_path.exists():
                zip_path.unlink()
            if extract_dir.exists():
                shutil.rmtree(extract_dir, ignore_errors=True)

    return DownloadArtifacts(
        genome_files=genome_files,
        warnings=warnings,
        notes=notes,
    )
