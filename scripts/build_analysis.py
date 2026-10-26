#!/usr/bin/env python3
"""Validate the transcribed DNS assay data and generate figures and report."""

from __future__ import annotations

import csv
import math
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib
from matplotlib.patches import FancyBboxPatch

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"
REPORTS = ROOT / "reports"

INK = "#17212B"
NAVY = "#17324D"
BLUE = "#3777A8"
ORANGE = "#D47A37"
TEAL = "#2A8C82"
MUTED = "#5F6B73"
GRID = "#DDE4E8"
PALE = "#EEF3F6"

ENZYME_COLORS = {
    "Commercial fungal amylase": BLUE,
    "Salivary amylase": ORANGE,
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def fit_standard_curve(points: list[dict[str, str]]) -> dict[str, float]:
    x = [float(row["maltose_concentration_mg_ml"]) for row in points]
    y = [float(row["absorbance_540nm"]) for row in points]
    x_mean = sum(x) / len(x)
    y_mean = sum(y) / len(y)
    denominator = sum((value - x_mean) ** 2 for value in x)
    if denominator == 0:
        raise ValueError("Standard concentrations must not all be identical.")
    slope = sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y)) / denominator
    intercept = y_mean - slope * x_mean
    residual_sum = sum((observed - (slope * concentration + intercept)) ** 2 for concentration, observed in zip(x, y))
    total_sum = sum((observed - y_mean) ** 2 for observed in y)
    r_squared = 1 - residual_sum / total_sum if total_sum else math.nan
    return {"slope": slope, "intercept": intercept, "r_squared": r_squared}


def read_project_data() -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, float]]:
    standards = read_csv(DATA / "maltose_standard_curve.csv")
    assays = read_csv(DATA / "assay_results.csv")
    if len(standards) != 9:
        raise ValueError(f"Expected 9 standard points; found {len(standards)}.")
    if len(assays) != 52:
        raise ValueError(f"Expected 52 assay readings; found {len(assays)}.")
    calibration = fit_standard_curve(standards)
    for row in assays:
        absorbance = float(row["absorbance_540nm"])
        if absorbance < 0:
            raise ValueError("Absorbance values must be nonnegative.")
    return standards, assays, calibration


def _set_style(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(GRID)
    ax.spines["bottom"].set_color(GRID)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, labelsize=9)


def _save_figure(fig: plt.Figure, stem: str, *, png_dpi: int = 300) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / f"{stem}.png", dpi=png_dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_standard_curve(
    standards: list[dict[str, str]], calibration: dict[str, float]
) -> None:
    x = [float(row["maltose_concentration_mg_ml"]) for row in standards]
    y = [float(row["absorbance_540nm"]) for row in standards]
    slope, intercept, r_squared = (
        calibration["slope"],
        calibration["intercept"],
        calibration["r_squared"],
    )
    line_x = [0, max(x)]
    line_y = [slope * value + intercept for value in line_x]

    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    fig.subplots_adjust(left=0.13, right=0.96, top=0.83, bottom=0.20)
    fig.suptitle(
        "Maltose calibration curve · DNS assay",
        x=0.13,
        y=0.96,
        ha="left",
        fontsize=15,
        fontweight="bold",
        color=INK,
    )
    ax.scatter(x, y, s=50, color=TEAL, edgecolor="white", linewidth=0.8, zorder=3, label="Reported standards")
    ax.plot(line_x, line_y, color=NAVY, linewidth=2, label="Ordinary least-squares fit")
    ax.set_xlabel("Maltose concentration (mg/mL)", color=MUTED)
    ax.set_ylabel("Absorbance at 540 nm", color=MUTED)
    _set_style(ax)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    sign = "−" if intercept < 0 else "+"
    equation = f"A540 = {slope:.5f} × C {sign} {abs(intercept):.5f}\nR² = {r_squared:.5f}"
    ax.text(
        0.98,
        0.06,
        equation,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color=INK,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": PALE, "edgecolor": "none"},
    )
    fig.text(
        0.13,
        0.07,
        "Fit recalculated from the nine tabulated standards; the source text contains an intercept-sign inconsistency.",
        fontsize=7.5,
        color=MUTED,
    )
    _save_figure(fig, "figure_01_maltose_standard_curve")


def plot_response(
    assays: list[dict[str, str]],
    *,
    experiment: str,
    title: str,
    xlabel: str,
    stem: str,
) -> None:
    selected = [row for row in assays if row["experiment"] == experiment]
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in selected:
        grouped[row["enzyme"]].append(row)

    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    fig.subplots_adjust(left=0.13, right=0.96, top=0.83, bottom=0.22)
    fig.suptitle(title, x=0.13, y=0.96, ha="left", fontsize=15, fontweight="bold", color=INK)
    tick_values = sorted(
        {float(row["condition_value"]) for row in selected}
    )
    for enzyme in ENZYME_COLORS:
        rows = sorted(grouped[enzyme], key=lambda row: float(row["condition_value"]))
        x = [float(row["condition_value"]) for row in rows]
        y = [float(row["absorbance_540nm"]) for row in rows]
        ax.plot(
            x,
            y,
            marker="o",
            markersize=6,
            linewidth=2,
            color=ENZYME_COLORS[enzyme],
            label=enzyme,
        )
    ax.set_xticks(tick_values)
    ax.set_xlabel(xlabel, color=MUTED)
    ax.set_ylabel("Reported absorbance at 540 nm", color=MUTED)
    _set_style(ax)
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.text(
        0.13,
        0.07,
        "One source-reported reading per condition; connecting lines are visual guides, not fitted models.",
        fontsize=7.5,
        color=MUTED,
    )
    _save_figure(fig, stem)


def plot_graphical_abstract(
    assays: list[dict[str, str]], calibration: dict[str, float]
) -> None:
    """Create a compact, source-grounded visual summary of all five tests."""
    specs = [
        ("pH", "pH", "pH", "peak"),
        ("temperature", "Temperature", "Temperature (°C)", "peak"),
        ("calcium_chloride", "CaCl₂", "CaCl₂ (source-labeled mM)", "peak"),
        ("starch_concentration", "Starch", "Starch (%; source-labeled)", "rise"),
        ("ethanol", "Ethanol", "Ethanol (%; final level unclear)", "fall"),
    ]
    fig = plt.figure(figsize=(15, 9.4), facecolor="#F5F8FA")
    fig.patches.append(
        plt.Rectangle(
            (0, 0.78),
            1,
            0.22,
            transform=fig.transFigure,
            facecolor=NAVY,
            edgecolor="none",
            zorder=-10,
        )
    )
    fig.text(0.055, 0.955, "GRAPHICAL ABSTRACT", fontsize=9, color="#A9D8D2", fontweight="bold")
    fig.text(
        0.055,
        0.895,
        "Amylase signal across physicochemical conditions",
        fontsize=22,
        color="white",
        fontweight="bold",
    )
    fig.text(
        0.055,
        0.825,
        "DNS endpoint A540 · Commercial fungal and salivary amylase",
        fontsize=10,
        color="#E7EEF2",
    )

    calibration_text = (
        "MALTOSE CALIBRATION REFIT\n"
        f"A540 = {calibration['slope']:.5f}C "
        f"{'−' if calibration['intercept'] < 0 else '+'} "
        f"{abs(calibration['intercept']):.5f}\n"
        f"R² = {calibration['r_squared']:.5f}"
    )
    fig.add_artist(
        FancyBboxPatch(
            (0.73, 0.815),
            0.22,
            0.14,
            boxstyle="round,pad=0.012,rounding_size=0.012",
            transform=fig.transFigure,
            facecolor="#264966",
            edgecolor="#486780",
            linewidth=0.8,
        )
    )
    fig.text(0.745, 0.94, calibration_text, fontsize=8, color="white", va="top", linespacing=1.35)

    legend_y = 0.752
    for x, color, label in [
        (0.058, ENZYME_COLORS["Commercial fungal amylase"], "Commercial fungal amylase"),
        (0.285, ENZYME_COLORS["Salivary amylase"], "Salivary amylase"),
    ]:
        fig.add_artist(
            plt.Line2D(
                [x, x + 0.025],
                [legend_y, legend_y],
                transform=fig.transFigure,
                color=color,
                marker="o",
                linewidth=2,
                markersize=4,
            )
        )
        fig.text(x + 0.032, legend_y, label, fontsize=8, color=INK, va="center")
    fig.text(0.945, legend_y, "Shared scale: 0–1.1 A540", fontsize=8, color=MUTED, ha="right", va="center")

    top_positions = [0.055, 0.365, 0.675]
    bottom_positions = [0.21, 0.52]
    card_width, card_height = 0.27, 0.295
    positions = [
        (top_positions[0], 0.435),
        (top_positions[1], 0.435),
        (top_positions[2], 0.435),
        (bottom_positions[0], 0.115),
        (bottom_positions[1], 0.115),
    ]

    for (experiment, title, x_label, trend), (card_x, card_y) in zip(specs, positions):
        selected = [row for row in assays if row["experiment"] == experiment]
        grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in selected:
            grouped[row["enzyme"]].append(row)
        ordered = {
            enzyme: sorted(rows, key=lambda row: float(row["condition_value"]))
            for enzyme, rows in grouped.items()
        }

        fig.add_artist(
            FancyBboxPatch(
                (card_x, card_y),
                card_width,
                card_height,
                boxstyle="round,pad=0.008,rounding_size=0.012",
                transform=fig.transFigure,
                facecolor="white",
                edgecolor=GRID,
                linewidth=0.9,
                zorder=-1,
            )
        )
        fig.text(
            card_x + 0.018,
            card_y + card_height - 0.035,
            title,
            fontsize=11,
            color=NAVY,
            fontweight="bold",
            va="center",
        )
        fig.text(
            card_x + card_width - 0.018,
            card_y + card_height - 0.035,
            x_label,
            fontsize=6.7,
            color=MUTED,
            ha="right",
            va="center",
        )

        ax = fig.add_axes([card_x + 0.045, card_y + 0.075, card_width - 0.07, card_height - 0.135])
        all_x = sorted({float(row["condition_value"]) for row in selected})
        for enzyme, rows in ordered.items():
            x = [float(row["condition_value"]) for row in rows]
            y = [float(row["absorbance_540nm"]) for row in rows]
            ax.plot(
                x,
                y,
                color=ENZYME_COLORS[enzyme],
                marker="o",
                markersize=3.5,
                linewidth=1.6,
                zorder=3,
            )
        ax.set_xticks(all_x)
        ax.set_xticklabels([f"{value:g}" for value in all_x])
        ax.set_ylim(0, 1.1)
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_ylabel("A540", fontsize=6.5, color=MUTED, labelpad=1)
        ax.tick_params(axis="both", colors=MUTED, labelsize=6.3, length=2, pad=1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color(GRID)
        ax.spines["bottom"].set_color(GRID)
        ax.grid(axis="y", color=GRID, linewidth=0.55)
        ax.set_axisbelow(True)
        ax.margins(x=0.04)

        if trend == "peak":
            peak_rows = [
                max(ordered[enzyme], key=lambda row: float(row["absorbance_540nm"]))
                for enzyme in ENZYME_COLORS
            ]
            condition = f"{float(peak_rows[0]['condition_value']):g}"
            summary = (
                f"Highest at {condition}: fungal {float(peak_rows[0]['absorbance_540nm']):.2f}; "
                f"salivary {float(peak_rows[1]['absorbance_540nm']):.2f}"
            )
        elif trend == "rise":
            fungal_end = float(ordered["Commercial fungal amylase"][-1]["absorbance_540nm"])
            saliva_end = float(ordered["Salivary amylase"][-1]["absorbance_540nm"])
            condition = f"{float(ordered['Commercial fungal amylase'][-1]['condition_value']):g}%"
            summary = f"Rises to {fungal_end:.2f} / {saliva_end:.2f} at {condition}; no plateau"
        else:
            fungal_end = float(ordered["Commercial fungal amylase"][-1]["absorbance_540nm"])
            saliva_end = float(ordered["Salivary amylase"][-1]["absorbance_540nm"])
            condition = f"{float(ordered['Commercial fungal amylase'][-1]['condition_value']):g}%"
            summary = f"Falls to {fungal_end:.2f} / {saliva_end:.2f} by {condition}"
        fig.text(
            card_x + 0.018,
            card_y + 0.027,
            textwrap.fill(summary, width=49),
            fontsize=7.2,
            color=INK,
            va="center",
        )

    fig.text(
        0.055,
        0.067,
        "One source-reported reading per condition; no replicates or error bars. Lines are visual guides, not fitted models.",
        fontsize=8.2,
        color=INK,
    )
    fig.text(
        0.055,
        0.038,
        "CaCl₂ and ethanol final concentrations are unclear. DNS results are endpoint absorbance, not normalized enzyme activity.",
        fontsize=7.6,
        color=MUTED,
    )
    fig.text(0.945, 0.038, "Source report: calibration pp. 8–9; assay tests pp. 13–32", fontsize=7, color=MUTED, ha="right")
    _save_figure(fig, "graphical_abstract", png_dpi=180)


def write_recalculated_yields(
    assays: list[dict[str, str]], calibration: dict[str, float]
) -> Path:
    RESULTS.mkdir(parents=True, exist_ok=True)
    output = RESULTS / "recalculated_product_yield.csv"
    slope, intercept = calibration["slope"], calibration["intercept"]
    columns = [
        "experiment",
        "condition_value",
        "condition_unit",
        "enzyme",
        "absorbance_540nm",
        "source_reported_maltose_mg_ml",
        "maltose_equivalent_from_refit_mg_ml",
        "source_minus_refit_mg_ml",
        "source_pages",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in assays:
            absorbance = float(row["absorbance_540nm"])
            concentration = (absorbance - intercept) / slope
            reported = float(row["source_reported_maltose_mg_ml"])
            writer.writerow(
                {
                    **{key: row[key] for key in columns[:6]},
                    "maltose_equivalent_from_refit_mg_ml": f"{concentration:.4f}",
                    "source_minus_refit_mg_ml": f"{reported - concentration:.4f}",
                    "source_pages": row["source_pages"],
                }
            )
    return output


def _text_page(pdf: PdfPages, number: int, title: str, sections: list[tuple[str, str]]) -> None:
    fig = plt.figure(figsize=(8.27, 11.69), facecolor="white")
    fig.text(0.08, 0.94, title, fontsize=20, fontweight="bold", color=INK, va="top")
    y = 0.88
    for heading, paragraph in sections:
        fig.text(0.08, y, heading, fontsize=12, fontweight="bold", color=NAVY, va="top")
        y -= 0.035
        for line in textwrap.fill(paragraph, width=96).splitlines():
            fig.text(0.08, y, line, fontsize=9, color=INK, va="top")
            y -= 0.021
        y -= 0.021
    _footer(fig, number)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _footer(fig: plt.Figure, page_number: int) -> None:
    fig.text(0.08, 0.035, "Amylase · DNS endpoint absorbance · Source-reported measurements", fontsize=7, color=MUTED)
    fig.text(0.92, 0.035, str(page_number), ha="right", fontsize=8, color=MUTED)


def _image_page(pdf: PdfPages, number: int, title: str, images: list[tuple[str, list[float]]]) -> None:
    fig = plt.figure(figsize=(8.27, 11.69), facecolor="white")
    fig.text(0.08, 0.94, title, fontsize=20, fontweight="bold", color=INK, va="top")
    for filename, box in images:
        ax = fig.add_axes(box)
        ax.imshow(plt.imread(FIGURES / filename))
        ax.axis("off")
    _footer(fig, number)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _graphical_abstract_page(pdf: PdfPages, number: int) -> None:
    fig = plt.figure(figsize=(11.69, 8.27), facecolor="white")
    ax = fig.add_axes([0.025, 0.075, 0.95, 0.87])
    ax.imshow(plt.imread(FIGURES / "graphical_abstract.png"))
    ax.axis("off")
    _footer(fig, number)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def build_pdf_report(calibration: dict[str, float]) -> Path:
    report_path = REPORTS / "final_report.pdf"
    REPORTS.mkdir(parents=True, exist_ok=True)
    with PdfPages(report_path) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69), facecolor="white")
        fig.text(0.08, 0.91, "ENZYMOLOGY · FINAL PROJECT REPORT", fontsize=10, color=TEAL, fontweight="bold")
        fig.text(
            0.08,
            0.82,
            "Characterization of amylase activity\nunder physicochemical conditions",
            fontsize=25,
            color=INK,
            fontweight="bold",
            va="top",
            linespacing=1.15,
        )
        fig.text(
            0.08,
            0.70,
            "DNS colorimetric assay · Commercial fungal and salivary amylase",
            fontsize=11,
            color=MUTED,
        )
        cards = [
            ("pH 7", "highest reported signal for both preparations"),
            ("37 °C", "highest reported signal in the temperature series"),
            ("10 mM", "source-labeled CaCl₂ response maximum"),
        ]
        for x, (value, caption) in zip([0.08, 0.39, 0.70], cards):
            fig.text(x, 0.59, value, fontsize=20, fontweight="bold", color=NAVY)
            fig.text(x, 0.54, textwrap.fill(caption, width=29), fontsize=8, color=MUTED, va="top")
            fig.add_artist(plt.Line2D([x, x + 0.23], [0.51, 0.51], transform=fig.transFigure, color=ORANGE, linewidth=3))
        fig.text(0.08, 0.42, "Executive summary", fontsize=14, fontweight="bold", color=NAVY)
        summary = (
            "The supplied report contains one A540 reading per condition across pH, temperature, "
            "CaCl₂, starch concentration, and ethanol. Signal rises through the tested starch range "
            "and falls as source-labeled ethanol increases. These are descriptive endpoint patterns, "
            "not replicated or normalized enzyme-rate estimates."
        )
        fig.text(0.08, 0.38, textwrap.fill(summary, width=96), fontsize=10, color=INK, va="top", linespacing=1.5)
        fig.text(0.08, 0.23, "Calibration refit", fontsize=12, fontweight="bold", color=NAVY)
        equation = (
            f"A540 = {calibration['slope']:.5f} × C "
            f"{'−' if calibration['intercept'] < 0 else '+'} "
            f"{abs(calibration['intercept']):.5f}; R² = {calibration['r_squared']:.5f}"
        )
        fig.text(0.08, 0.19, equation, fontsize=10, color=INK)
        fig.text(
            0.08,
            0.12,
            "No replicate readings or raw instrument export were supplied. Interpret results as preliminary descriptive evidence.",
            fontsize=8,
            color=MUTED,
        )
        _footer(fig, 1)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        _graphical_abstract_page(pdf, 2)
        _image_page(
            pdf,
            3,
            "Maltose calibration",
            [("figure_01_maltose_standard_curve.png", [0.09, 0.17, 0.82, 0.68])],
        )
        _image_page(
            pdf,
            4,
            "pH and temperature responses",
            [
                ("figure_02_ph_response.png", [0.07, 0.52, 0.86, 0.31]),
                ("figure_03_temperature_response.png", [0.07, 0.13, 0.86, 0.31]),
            ],
        )
        _image_page(
            pdf,
            5,
            "Calcium and substrate responses",
            [
                ("figure_04_calcium_chloride_response.png", [0.07, 0.52, 0.86, 0.31]),
                ("figure_05_starch_response.png", [0.07, 0.13, 0.86, 0.31]),
            ],
        )
        _image_page(
            pdf,
            6,
            "Ethanol response",
            [("figure_06_ethanol_response.png", [0.08, 0.23, 0.84, 0.60])],
        )
        _text_page(
            pdf,
            7,
            "Interpretation and next steps",
            [
                (
                    "What the measurements support",
                    "Within the supplied one-reading-per-condition data, both preparations show their "
                    "highest pH signal at 7 and highest temperature signal at 37°C. Both peak at the "
                    "source-labeled 10 mM CaCl₂ condition. The starch response increases through 0.5%, "
                    "while ethanol-associated absorbance decreases through 1.5%.",
                ),
                (
                    "What they do not establish",
                    "There are no replicate values or uncertainty estimates, so the report does not "
                    "support significance testing. The starch series does not demonstrate a plateau; "
                    "Km and Vmax cannot be inferred. Different enzyme preparations are not normalized "
                    "to protein or activity units.",
                ),
                (
                    "Assay and reporting caveats",
                    "DNS responds to reducing sugars, so calibration output is maltose-equivalent. "
                    "Buffer chemistry changes across the pH series; final CaCl₂ and ethanol concentrations "
                    "are ambiguous. The source calibration text contains opposite intercept signs, and "
                    "some ethanol concentrations do not match the refitted equation within rounding.",
                ),
                (
                    "Recommended follow-up",
                    "Repeat with independent replicates, matched blanks and controls, recorded final "
                    "reaction concentrations, and normalized enzyme input. Preserve raw spectrophotometer "
                    "exports and refit the standard curve from replicate standards.",
                ),
            ],
        )
    return report_path


def build_analysis() -> None:
    standards, assays, calibration = read_project_data()
    plot_standard_curve(standards, calibration)
    response_specs = [
        ("pH", "Effect of pH on amylase signal", "pH", "figure_02_ph_response"),
        ("temperature", "Temperature response of amylase", "Temperature (°C)", "figure_03_temperature_response"),
        (
            "calcium_chloride",
            "CaCl₂ response of amylase",
            "Source-labeled CaCl₂ concentration (mM; final level unclear)",
            "figure_04_calcium_chloride_response",
        ),
        (
            "starch_concentration",
            "Starch concentration response",
            "Source-labeled starch concentration (%)",
            "figure_05_starch_response",
        ),
        (
            "ethanol",
            "Ethanol response of amylase",
            "Source-labeled ethanol concentration (%; final level unclear)",
            "figure_06_ethanol_response",
        ),
    ]
    for experiment, title, xlabel, stem in response_specs:
        plot_response(assays, experiment=experiment, title=title, xlabel=xlabel, stem=stem)
    plot_graphical_abstract(assays, calibration)
    write_recalculated_yields(assays, calibration)
    report_path = build_pdf_report(calibration)
    print(f"Calibration: A540 = {calibration['slope']:.5f} C {calibration['intercept']:+.5f}")
    print(f"R^2 = {calibration['r_squared']:.5f}; assay readings = {len(assays)}")
    print(f"Generated figures in {FIGURES}")
    print(f"Generated report: {report_path}")


if __name__ == "__main__":
    build_analysis()
