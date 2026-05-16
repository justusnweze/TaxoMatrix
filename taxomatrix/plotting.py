from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from matplotlib import cm
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform


NEUTRAL_METADATA_COLOR = "#d9d9d9"


@dataclass(slots=True)
class HeatmapArtifacts:
    output_paths: list[Path]
    metadata_legend_path: Path | None
    metadata_sheet: pd.DataFrame | None
    warnings: list[str]


def _build_distance_matrix(ani_matrix: pd.DataFrame) -> pd.DataFrame:
    filled = ani_matrix.fillna(0.0).clip(lower=0.0, upper=100.0)
    distance = 100.0 - filled
    for label in distance.index:
        distance.loc[label, label] = 0.0
    return distance


def _build_linkage(ani_matrix: pd.DataFrame):
    distance = _build_distance_matrix(ani_matrix)
    condensed = squareform(distance.to_numpy(), checks=False)
    return linkage(condensed, method="average")


def _determine_rotation(genome_count: int) -> int:
    if genome_count <= 6:
        return 40
    if genome_count <= 20:
        return 60
    return 90


def _determine_font_size(genome_count: int, max_label_length: int) -> float:
    if genome_count <= 6 and max_label_length <= 20:
        return 10.0
    if genome_count <= 10:
        return 9.0
    if genome_count <= 20:
        return 8.0
    if genome_count <= 40:
        return 7.0
    return 6.0


def _compute_figure_dimensions(genome_count: int, max_label_length: int, metadata_tracks: int) -> tuple[float, float]:
    base_heatmap = min(max(4.3, genome_count * 0.45), 14.5)
    width = base_heatmap + min(max(2.4, max_label_length * 0.075), 6.5) + (metadata_tracks * 0.25) + 1.7
    height = base_heatmap + min(max(1.8, max_label_length * 0.04), 4.8) + 1.4
    return width, height


def _should_show_values(show_values: bool) -> bool:
    return show_values


def _coerce_metadata_value(value) -> str:
    if pd.isna(value):
        return "Missing"
    text = str(value).strip()
    return text if text else "Missing"


def _build_categorical_color_map(series: pd.Series, palette_name: str) -> dict[str, str]:
    values = sorted({_coerce_metadata_value(value) for value in series})
    palette = sns.color_palette(palette_name, n_colors=max(1, len(values)))
    return {value: mcolors.to_hex(color) for value, color in zip(values, palette, strict=True)}


def _load_metadata(
    metadata_path: Path,
    label_mapping: pd.DataFrame,
    reordered_display_labels: list[str],
    reordered_figure_labels: list[str],
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, Path, pd.DataFrame, list[str]]:
    metadata = pd.read_csv(metadata_path, sep="\t")
    required_columns = {"genome", "habitat", "source", "pH"}
    missing_columns = sorted(required_columns - set(metadata.columns))
    if missing_columns:
        raise ValueError(
            "Metadata file must contain these columns: genome, habitat, source, pH. "
            f"Missing: {', '.join(missing_columns)}"
        )

    label_lookup = dict(zip(label_mapping["display_label"], label_mapping["display_label"], strict=True))
    accession_lookup = dict(zip(label_mapping["assembly_accession"], label_mapping["display_label"], strict=True))
    figure_lookup = dict(zip(label_mapping["display_label"], reordered_figure_labels, strict=True))

    metadata = metadata.copy()
    metadata["matched_display_label"] = metadata["genome"].map(label_lookup)
    unresolved = metadata["matched_display_label"].isna()
    metadata.loc[unresolved, "matched_display_label"] = metadata.loc[unresolved, "genome"].map(accession_lookup)

    unmatched_values = metadata.loc[metadata["matched_display_label"].isna(), "genome"].tolist()
    warnings: list[str] = []
    if unmatched_values:
        warnings.append(
            "Metadata rows could not be matched to genomes and were ignored: "
            + ", ".join(map(str, unmatched_values))
        )

    metadata = metadata.dropna(subset=["matched_display_label"]).drop_duplicates(
        subset=["matched_display_label"],
        keep="first",
    )
    metadata = metadata.set_index("matched_display_label")
    metadata = metadata.reindex(reordered_display_labels)
    metadata["figure_label"] = metadata.index.map(figure_lookup)

    legend_records: list[dict[str, str]] = []
    row_color_columns: list[pd.Series] = []

    for field, palette_name in (("habitat", "Set2"), ("source", "tab10")):
        color_map = _build_categorical_color_map(metadata[field], palette_name)
        color_map["Missing"] = NEUTRAL_METADATA_COLOR
        row_color_columns.append(
            metadata[field].map(_coerce_metadata_value).map(color_map).set_axis(reordered_figure_labels).rename(field)
        )
        missing_count = int(metadata[field].isna().sum())
        if missing_count:
            warnings.append(
                f"Metadata field '{field}' was missing for {missing_count} genome(s); neutral colors were used."
            )
        for value, color in color_map.items():
            legend_records.append({"field": field, "value": value, "color": color})

    ph_series = pd.to_numeric(metadata["pH"], errors="coerce")
    ph_missing = int(ph_series.isna().sum())
    if ph_missing:
        warnings.append(f"Metadata field 'pH' was missing for {ph_missing} genome(s); neutral colors were used.")

    valid_ph = ph_series.dropna()
    cmap = cm.get_cmap("YlGnBu")
    if valid_ph.empty:
        ph_colors = [NEUTRAL_METADATA_COLOR] * len(metadata.index)
    else:
        if float(valid_ph.min()) == float(valid_ph.max()):
            ph_colors = [mcolors.to_hex(cmap(0.5)) if pd.notna(value) else NEUTRAL_METADATA_COLOR for value in ph_series]
        else:
            norm = mcolors.Normalize(vmin=float(valid_ph.min()), vmax=float(valid_ph.max()))
            ph_colors = [
                mcolors.to_hex(cmap(norm(float(value)))) if pd.notna(value) else NEUTRAL_METADATA_COLOR
                for value in ph_series
            ]
        legend_records.extend(
            [
                {"field": "pH", "value": f"min={float(valid_ph.min()):.2f}", "color": mcolors.to_hex(cmap(0.0))},
                {"field": "pH", "value": f"max={float(valid_ph.max()):.2f}", "color": mcolors.to_hex(cmap(1.0))},
                {"field": "pH", "value": "Missing", "color": NEUTRAL_METADATA_COLOR},
            ]
        )
    row_color_columns.append(pd.Series(ph_colors, index=reordered_figure_labels, name="pH"))

    row_colors = pd.concat(row_color_columns, axis=1)
    legend_frame = pd.DataFrame.from_records(legend_records)
    legend_path = output_dir / "metadata_legend.tsv"
    legend_frame.to_csv(legend_path, sep="\t", index=False)

    metadata_sheet = metadata.reset_index(names="display_label")
    ordered_columns = ["figure_label", "display_label", "genome", "habitat", "source", "pH"]
    metadata_sheet = metadata_sheet.loc[:, [column for column in ordered_columns if column in metadata_sheet.columns]]

    return metadata, row_colors, legend_path, metadata_sheet, warnings


def _set_annotation_contrast(cluster_grid, data: pd.DataFrame, cmap_name: str) -> None:
    cmap = cm.get_cmap(cmap_name)
    norm = mcolors.Normalize(vmin=0.0, vmax=100.0)
    texts = cluster_grid.ax_heatmap.texts
    values = data.to_numpy().flatten()
    for text, value in zip(texts, values, strict=False):
        r, g, b, _ = cmap(norm(float(value)))
        luminance = (0.2126 * r) + (0.7152 * g) + (0.0722 * b)
        text.set_color("black" if luminance > 0.62 else "white")


def plot_clustered_heatmap(
    ani_matrix: pd.DataFrame,
    output_dir: Path,
    *,
    label_mapping: pd.DataFrame,
    show_values: bool = False,
    species_boundary: float = 95.0,
    annotate_boundary: bool = False,
    metadata_path: Path | None = None,
) -> HeatmapArtifacts:
    if ani_matrix.shape[0] < 2:
        raise ValueError("At least two genomes are required to plot a clustered heatmap.")

    sns.set_theme(style="white")
    output_prefix = output_dir / "ani_clustered_heatmap"
    warnings: list[str] = []

    figure_lookup = dict(zip(label_mapping["display_label"], label_mapping["short_figure_label"], strict=True))
    reordered_display_labels = ani_matrix.index.tolist()
    figure_labels = [figure_lookup[label] for label in reordered_display_labels]

    figure_matrix = ani_matrix.rename(index=figure_lookup, columns=figure_lookup)
    plot_data = figure_matrix.fillna(0.0)
    mask = figure_matrix.isna()

    genome_count = figure_matrix.shape[0]
    max_label_length = max(len(label) for label in figure_matrix.index)
    row_linkage = _build_linkage(figure_matrix)
    col_linkage = _build_linkage(figure_matrix.T)

    row_colors = None
    metadata_legend_path = None
    metadata_sheet = None
    if metadata_path is not None:
        _, row_colors, metadata_legend_path, metadata_sheet, metadata_warnings = _load_metadata(
            metadata_path=metadata_path,
            label_mapping=label_mapping,
            reordered_display_labels=reordered_display_labels,
            reordered_figure_labels=figure_labels,
            output_dir=output_dir,
        )
        warnings.extend(metadata_warnings)

    rotation = _determine_rotation(genome_count)
    font_size = _determine_font_size(genome_count, max_label_length)
    figure_width, figure_height = _compute_figure_dimensions(
        genome_count=genome_count,
        max_label_length=max_label_length,
        metadata_tracks=0 if row_colors is None else row_colors.shape[1],
    )
    annot = _should_show_values(show_values)
    annot_font_size = max(4.0, font_size - 1.5)

    cluster_grid = sns.clustermap(
        plot_data,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        row_colors=row_colors,
        col_colors=row_colors,
        cmap="viridis",
        vmin=0,
        vmax=100,
        linewidths=0.35,
        linecolor="white",
        mask=mask,
        figsize=(figure_width, figure_height),
        dendrogram_ratio=(0.10, 0.10),
        colors_ratio=(0.025 if row_colors is not None else 0.01, 0.025 if row_colors is not None else 0.01),
        cbar_pos=(0.82, 0.74, 0.022, 0.18),
        cbar_kws={"label": "ANI (%)"},
        xticklabels=True,
        yticklabels=True,
        annot=annot,
        fmt=".1f",
        annot_kws={"fontsize": annot_font_size},
    )

    cluster_grid.fig.subplots_adjust(left=0.04, right=0.94, bottom=0.08, top=0.96)
    cluster_grid.ax_heatmap.set_position([0.20, 0.14, 0.62, 0.68])
    cluster_grid.ax_row_dendrogram.set_position([0.08, 0.14, 0.10, 0.68])
    cluster_grid.ax_col_dendrogram.set_position([0.20, 0.84, 0.62, 0.10])
    if cluster_grid.ax_row_colors is not None:
        cluster_grid.ax_row_colors.set_position([0.18, 0.14, 0.016, 0.68])
        cluster_grid.ax_row_colors.set_xticks([])
        cluster_grid.ax_row_colors.set_yticks([])
    if cluster_grid.ax_col_colors is not None:
        cluster_grid.ax_col_colors.set_position([0.20, 0.815, 0.62, 0.018])
        cluster_grid.ax_col_colors.set_xticks([])
        cluster_grid.ax_col_colors.set_yticks([])
    if cluster_grid.cax is not None:
        cluster_grid.cax.set_position([0.84, 0.71, 0.024, 0.19])
        cluster_grid.cax.tick_params(labelsize=max(6.0, font_size - 1))
        cluster_grid.cax.yaxis.label.set_size(font_size)

    cluster_grid.ax_heatmap.set_xlabel("Genome", fontsize=font_size + 1)
    cluster_grid.ax_heatmap.set_ylabel("Genome", fontsize=font_size + 1)
    cluster_grid.ax_heatmap.tick_params(axis="x", labelsize=font_size, pad=1)
    cluster_grid.ax_heatmap.tick_params(axis="y", labelsize=font_size, pad=1)
    plt.setp(cluster_grid.ax_heatmap.get_xticklabels(), rotation=rotation, ha="right", rotation_mode="anchor")
    plt.setp(cluster_grid.ax_heatmap.get_yticklabels(), rotation=0)

    if annot:
        _set_annotation_contrast(cluster_grid, plot_data, "viridis")

    if annotate_boundary and cluster_grid.cax is not None:
        cluster_grid.fig.text(
            0.80,
            0.67,
            f"Reference species boundary: {species_boundary:.1f}% ANI",
            ha="left",
            va="top",
            fontsize=max(6.0, font_size - 1),
        )
        cluster_grid.cax.axhline(species_boundary, color="black", linestyle="--", linewidth=0.9, alpha=0.8)

    output_paths: list[Path] = []
    for suffix in ("png", "svg", "pdf"):
        path = output_prefix.with_suffix(f".{suffix}")
        save_kwargs = {"bbox_inches": "tight", "pad_inches": 0.08}
        if suffix == "png":
            save_kwargs["dpi"] = 300
        cluster_grid.fig.savefig(path, **save_kwargs)
        output_paths.append(path)

    plt.close(cluster_grid.fig)
    return HeatmapArtifacts(
        output_paths=output_paths,
        metadata_legend_path=metadata_legend_path,
        metadata_sheet=metadata_sheet,
        warnings=warnings,
    )
