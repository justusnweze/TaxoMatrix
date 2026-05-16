from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tempfile

from taxomatrix.utils import ensure_directory, resolve_executable, run_command


@dataclass(slots=True)
class FastANIConfig:
    executable: str = "fastANI"
    threads: int = 1
    frag_len: int | None = None
    min_fraction: float | None = None


@dataclass(slots=True)
class FastANIResult:
    raw_output_path: Path
    executable_path: str
    command: list[str]


def get_fastani_version(executable: str) -> str:
    executable_path = resolve_executable(executable)
    completed = run_command([executable_path, "--version"])
    version = completed.stdout.strip() or completed.stderr.strip()
    return version or "unknown"


def run_all_vs_all(
    genome_files: list[Path],
    output_dir: Path,
    config: FastANIConfig,
) -> FastANIResult:
    if len(genome_files) < 2:
        raise ValueError(
            "Fewer than 2 genome FASTA files were provided. "
            "At least 2 genomes are required for all-vs-all ANI."
        )

    ensure_directory(output_dir)
    executable = resolve_executable(config.executable)
    raw_output_path = output_dir / "fastani_raw.tsv"

    with tempfile.TemporaryDirectory(prefix="taxomatrix_fastani_", dir=output_dir) as temp_dir:
        temp_dir_path = Path(temp_dir)
        query_list = temp_dir_path / "queries.txt"
        ref_list = temp_dir_path / "references.txt"

        query_list.write_text(
            "\n".join(str(path.resolve()) for path in genome_files) + "\n",
            encoding="utf-8",
        )
        ref_list.write_text(
            "\n".join(str(path.resolve()) for path in genome_files) + "\n",
            encoding="utf-8",
        )

        command = [
            executable,
            "--ql",
            str(query_list),
            "--rl",
            str(ref_list),
            "-o",
            str(raw_output_path),
            "-t",
            str(config.threads),
        ]
        if config.frag_len is not None:
            command.extend(["--fragLen", str(config.frag_len)])
        if config.min_fraction is not None:
            command.extend(["--minFraction", str(config.min_fraction)])

        run_command(command, cwd=output_dir)

    if not raw_output_path.exists():
        raise FileNotFoundError(f"FastANI did not produce the expected output file: {raw_output_path}")

    return FastANIResult(
        raw_output_path=raw_output_path,
        executable_path=executable,
        command=command,
    )
