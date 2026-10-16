"""Reusable Matplotlib helpers for the scientific-line-plot house style."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib as mpl
from matplotlib import font_manager
from matplotlib.axes import Axes
from matplotlib.figure import Figure


@dataclass(frozen=True)
class ScientificLineStyle:
    axis_label_size: float = 22.0
    legend_size: float = 20.0
    tick_label_size: float = 20.0
    line_width: float = 2.0
    marker_size: float = 6.0
    spine_width: float = 1.1
    tick_width: float = 1.0
    grid_width: float = 0.7
    grid_color: str = "#B8C0C8"
    grid_alpha: float = 0.35
    png_dpi: int = 320


STYLE = ScientificLineStyle()


def _register_windows_times_new_roman() -> None:
    font_dir = Path(r"C:\Windows\Fonts")
    for filename in ("times.ttf", "timesbd.ttf", "timesi.ttf", "timesbi.ttf"):
        path = font_dir / filename
        if path.exists():
            font_manager.fontManager.addfont(path)


def configure_matplotlib(style: ScientificLineStyle = STYLE) -> str:
    """Configure Matplotlib and return the resolved Times New Roman family name.

    Fail explicitly when Times New Roman is unavailable instead of silently using
    a different publication font.
    """

    _register_windows_times_new_roman()
    family = "Times New Roman"
    try:
        font_manager.findfont(
            font_manager.FontProperties(family=family), fallback_to_default=False
        )
    except ValueError as exc:
        raise RuntimeError(
            "Times New Roman is not installed. Install or supply it before using "
            "this house style, or follow an explicit user-approved substitute."
        ) from exc

    mpl.rcParams.update(
        {
            "font.family": family,
            "font.weight": "bold",
            "axes.labelweight": "bold",
            "mathtext.fontset": "stix",
            "mathtext.default": "bf",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )
    return family


def style_axes(
    axis: Axes,
    *,
    xlabel: str,
    ylabel: str,
    style: ScientificLineStyle = STYLE,
    show_major_grid: bool = True,
) -> None:
    """Apply symbol-only labels, boxed spines, and bottom/left ticks."""

    axis.set_xlabel(xlabel, fontsize=style.axis_label_size, fontweight="bold")
    axis.set_ylabel(ylabel, fontsize=style.axis_label_size, fontweight="bold")
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(style.spine_width)

    axis.tick_params(
        axis="both",
        which="both",
        direction="out",
        bottom=True,
        left=True,
        top=False,
        right=False,
        labeltop=False,
        labelright=False,
        labelsize=style.tick_label_size,
        width=style.tick_width,
    )
    for label in list(axis.get_xticklabels()) + list(axis.get_yticklabels()):
        label.set_fontweight("bold")

    axis.grid(False, which="minor")
    axis.grid(
        show_major_grid,
        which="major",
        color=style.grid_color,
        linewidth=style.grid_width,
        alpha=style.grid_alpha,
    )


def style_legend(
    axis: Axes,
    *,
    handles: Sequence | None = None,
    labels: Sequence[str] | None = None,
    loc: str = "best",
    ncol: int = 1,
    style: ScientificLineStyle = STYLE,
    **kwargs,
):
    """Create a frameless bold Times New Roman legend."""

    prop = font_manager.FontProperties(
        family="Times New Roman", weight="bold", size=style.legend_size
    )
    legend_options = {
        "loc": loc,
        "ncol": ncol,
        "frameon": False,
        "prop": prop,
        **kwargs,
    }
    if handles is not None:
        legend_options["handles"] = handles
    if labels is not None:
        legend_options["labels"] = labels
    legend = axis.legend(**legend_options)
    for text in legend.get_texts():
        text.set_fontweight("bold")
    return legend


def style_data_lines(
    axis: Axes,
    *,
    style: ScientificLineStyle = STYLE,
    skip_labels: Iterable[str] = (),
) -> None:
    """Normalize existing data lines to 2 pt without changing colors or markers."""

    skipped = set(skip_labels)
    for line in axis.lines:
        if line.get_label() in skipped:
            continue
        line.set_linewidth(style.line_width)


def save_png_pdf(
    figure: Figure,
    output_stem: str | Path,
    *,
    style: ScientificLineStyle = STYLE,
) -> tuple[Path, Path]:
    """Save matching 320 dpi PNG and vector PDF using atomic replacement."""

    stem = Path(output_stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    temporary_png = stem.parent / f"{stem.name}.new.png"
    temporary_pdf = stem.parent / f"{stem.name}.new.pdf"
    figure.savefig(temporary_png, dpi=style.png_dpi)
    figure.savefig(temporary_pdf)
    temporary_png.replace(png)
    temporary_pdf.replace(pdf)
    return png, pdf
