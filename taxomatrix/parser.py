from __future__ import annotations

from pathlib import Path

import pandas as pd

FASTANI_COLUMNS = [
    "query",
    "reference",
    "ani",
    "fragments_matched",
    "total_fragments",
]


def parse_fastani_output(raw_output_path: Path) -> pd.DataFrame:
    if raw_output_path.stat().st_size == 0:
        raise ValueError("FastANI output is empty. Check your input genomes and FastANI run.")

    frame = pd.read_csv(
        raw_output_path,
        sep="\t",
        header=None,
        names=FASTANI_COLUMNS,
    )
    frame["query"] = frame["query"].map(lambda value: Path(value))
    frame["reference"] = frame["reference"].map(lambda value: Path(value))
    return frame


def build_ani_matrix(
    results: pd.DataFrame,
    genome_files: list[Path],
    path_to_label: dict[Path, str],
) -> pd.DataFrame:
    labels = [path_to_label[path.resolve()] for path in genome_files]

    matrix = pd.DataFrame(index=labels, columns=labels, dtype=float)
    for label in labels:
        matrix.loc[label, label] = 100.0

    for row in results.itertuples(index=False):
        query_label = path_to_label.get(Path(row.query).resolve())
        reference_label = path_to_label.get(Path(row.reference).resolve())
        if query_label is None or reference_label is None:
            continue

        current = matrix.loc[query_label, reference_label]
        new_value = float(row.ani)
        if pd.isna(current) or new_value > current:
            matrix.loc[query_label, reference_label] = new_value
            matrix.loc[reference_label, query_label] = new_value

    matrix = matrix.sort_index(axis=0).sort_index(axis=1)
    return matrix
