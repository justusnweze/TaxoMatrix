from __future__ import annotations

from datetime import datetime
from pathlib import Path
import re
import shlex
import shutil
import subprocess
from typing import Iterable

import pandas as pd


FASTA_SUFFIXES = {
    ".fa",
    ".fna",
    ".fasta",
    ".fas",
    ".fsa",
    ".ffn",
    ".faa",
}

ASSEMBLY_ACCESSION_PATTERN = re.compile(r"(GC[AF]_\d+\.\d+)")
HEADER_TAG_PATTERN = re.compile(r"\[(\w+)=([^\]]+)\]")
STRAIN_PATTERN = re.compile(r"\b(strain|isolate)\b\s+(.+)$", re.IGNORECASE)
SEQUENCE_ACCESSION_PATTERN = re.compile(r"^[A-Z]{1,4}_?[A-Z0-9]*\d[\w.:-]*$", re.IGNORECASE)
DISPLAY_LABEL_ACCESSION_PATTERN = re.compile(r"\s*\((GC[AF]_\d+\.\d+)\)\s*$")
TRAILING_METADATA_KEYWORDS = [
    " whole genome shotgun sequence",
    " complete genome",
    " complete sequence",
    " chromosome",
    " genomic scaffold",
    " scaffold",
    " contig",
]


def ensure_directory(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def discover_fasta_files(input_dir: Path) -> list[Path]:
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input path is not a directory: {input_dir}")

    fasta_files = sorted(
        path
        for path in input_dir.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
    )
    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA files were found in {input_dir}. "
            "TaxoMatrix requires genome FASTA files with one of these extensions: "
            f"{', '.join(sorted(FASTA_SUFFIXES))}"
        )
    return fasta_files


def genome_label(path: Path) -> str:
    return path.stem


def extract_assembly_accession(path: Path) -> str:
    match = ASSEMBLY_ACCESSION_PATTERN.search(path.name)
    if not match:
        raise ValueError(
            f"Could not determine a valid assembly accession from filename: {path.name}. "
            "Expected an accession such as GCA_000000000.1 or GCF_000000000.1 in the FASTA filename."
        )
    return match.group(1)


def read_first_fasta_header(path: Path) -> str | None:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                return line[1:].strip()
    return None


def _clean_header_title(header: str) -> str:
    stripped = HEADER_TAG_PATTERN.sub("", header).strip()
    if " " in stripped:
        first_token, remainder = stripped.split(" ", 1)
        if SEQUENCE_ACCESSION_PATTERN.match(first_token):
            stripped = remainder.strip()

    lowered = stripped.lower()
    cut_positions = [lowered.find(keyword) for keyword in TRAILING_METADATA_KEYWORDS if lowered.find(keyword) > 0]
    if cut_positions:
        stripped = stripped[: min(cut_positions)].strip(" ,;")

    return re.sub(r"\s+", " ", stripped).strip(" ,;")


def _build_display_core(organism_name: str, strain: str, cleaned_title: str) -> str:
    if cleaned_title:
        return cleaned_title
    if organism_name and strain:
        if strain.lower() in organism_name.lower():
            return organism_name
        return f"{organism_name} {strain}"
    return organism_name


def build_display_label(display_core: str, assembly_accession: str) -> str:
    return f"{display_core} ({assembly_accession})" if display_core else assembly_accession


def build_short_figure_label(display_label: str) -> str:
    accession_match = DISPLAY_LABEL_ACCESSION_PATTERN.search(display_label)
    core = DISPLAY_LABEL_ACCESSION_PATTERN.sub("", display_label).strip()
    if not core:
        return display_label

    parts = core.split()
    if len(parts) >= 2 and parts[0][0].isalpha():
        abbreviated_core = " ".join([f"{parts[0][0]}.", *parts[1:]])
        return abbreviated_core
    return core


def build_genome_label_mapping(genome_files: list[Path]) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    for genome_file in genome_files:
        assembly_accession = extract_assembly_accession(genome_file)
        header = read_first_fasta_header(genome_file)

        organism_name = ""
        strain = ""
        cleaned_title = ""

        if header:
            tags = {key.lower(): value.strip() for key, value in HEADER_TAG_PATTERN.findall(header)}
            organism_name = tags.get("organism", "").strip()
            strain = (tags.get("strain") or tags.get("isolate") or "").strip()
            cleaned_title = _clean_header_title(header)

            if cleaned_title:
                strain_match = STRAIN_PATTERN.search(cleaned_title)
                if strain_match:
                    parsed_organism = cleaned_title[: strain_match.start()].strip(" ,;")
                    parsed_strain = strain_match.group(2).strip(" ,;")
                    if not organism_name:
                        organism_name = parsed_organism
                    if not strain:
                        strain = parsed_strain
                elif not organism_name:
                    organism_name = cleaned_title

        display_core = _build_display_core(organism_name, strain, cleaned_title)
        display_label = build_display_label(display_core, assembly_accession)
        figure_label = build_short_figure_label(display_label)

        records.append(
            {
                "original_filename": genome_file.name,
                "assembly_accession": assembly_accession,
                "organism_name": organism_name,
                "strain": strain,
                "display_label": display_label,
                "short_figure_label": figure_label,
            }
        )

    return pd.DataFrame.from_records(records)


def validate_unique_genome_labels(genome_files: list[Path]) -> None:
    labels = [genome_label(path) for path in genome_files]
    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise ValueError(
            "Duplicate genome names detected after removing file extensions. "
            "Each genome filename stem must be unique. "
            f"Duplicates: {', '.join(duplicates)}"
        )


def validate_unique_display_labels(label_mapping: pd.DataFrame) -> None:
    labels = label_mapping["display_label"].tolist()
    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise ValueError(
            "Duplicate publication-friendly genome labels were generated. "
            "Check the metadata or filenames for collisions. "
            f"Duplicates: {', '.join(duplicates)}"
        )


def validate_unique_figure_labels(label_mapping: pd.DataFrame) -> None:
    labels = label_mapping["short_figure_label"].tolist()
    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise ValueError(
            "Duplicate short figure labels were generated. "
            "Adjust the underlying genome metadata so the abbreviated figure labels remain unique. "
            f"Duplicates: {', '.join(duplicates)}"
        )


def resolve_executable(executable: str) -> str:
    resolved = shutil.which(executable)
    if resolved:
        return resolved

    candidate = Path(executable)
    if candidate.exists() and candidate.is_file():
        return str(candidate)

    raise FileNotFoundError(
        f"FastANI executable not found: {executable}. "
        "Install FastANI or provide its full path with --fastani-binary."
    )


def run_command(command: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=str(cwd) if cwd else None,
        check=True,
        text=True,
        capture_output=True,
    )


def format_command(command: list[str]) -> str:
    return shlex.join(command)


def write_run_log(
    output_path: Path,
    *,
    command: str,
    fastani_command: str | None,
    runtime_seconds: float,
    genome_count: int,
    fastani_version: str,
    genome_labels: Iterable[str],
    output_paths: Iterable[str],
    warnings: Iterable[str],
) -> Path:
    warning_lines = list(warnings)
    label_lines = list(genome_labels)
    output_lines = list(output_paths)
    log_lines = [
        f"timestamp: {datetime.now().isoformat(timespec='seconds')}",
        f"command: {command}",
        f"runtime_seconds: {runtime_seconds:.2f}",
        f"genome_count: {genome_count}",
        f"fastani_version: {fastani_version}",
    ]
    if fastani_command:
        log_lines.append(f"fastani_command: {fastani_command}")
    log_lines.append("genome_labels:")
    log_lines.extend(f"- {label}" for label in label_lines)
    log_lines.append("output_paths:")
    log_lines.extend(f"- {path}" for path in output_lines)
    log_lines.append("warnings:")
    if warning_lines:
        log_lines.extend(f"- {warning}" for warning in warning_lines)
    else:
        log_lines.append("- none")

    output_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    return output_path
