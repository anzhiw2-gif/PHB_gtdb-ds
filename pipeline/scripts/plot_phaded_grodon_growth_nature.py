#!/usr/bin/env python3
"""Nature-style figure for the PhaDED gRodon2 growth-rate comparison.

Style follows the project's nature-figure convention
(``nature_figures_common.py`` / ``nature_figure_growth_comparison.py``):
Arial-first sans-serif, ``svg.fonttype="none"`` (editable SVG text),
``pdf.fonttype=42``, 7.2-inch double-column width, fixed palette, bold 9 pt
panel letters with an overlap QC pass, SVG/PDF/PNG(600 dpi) export, and
source-data tables written next to the figure.

Panels
  a  Genome-level predicted growth rate: candidate-carrier genomes vs their
     same-genus controls with no candidate detected under the defined search
     (box + jitter, Mann-Whitney U and Cliff's delta; descriptive only, and
     labelled secondary because genome-level sampling is not independent).
  b  Genus-level delta distribution (analysis unit = genus) with the
     signed-rank / stratified-permutation / sign-test summary.
  c  Paired genus means with the 1:1 line (colour = delta, size = pairs).
  d  Intracellular vs extracellular candidate carriers: genus-delta
     distributions and the between-group test.  **Secondary analysis**: the
     compartment comparison is not the primary estimand.
  e  Superfamily forest plot: mean genus delta with bootstrap 95% CI and
     Benjamini-Hochberg q values.
  f  Data ledger (accounting boxes) and the headline statements.

Group labels are "candidate gene carrier" and "candidate not detected under
the defined search"; the figure never emits a phenotype label.  The wording
for a null result is "no difference detected under the current design", and
equivalence wording appears only when an equivalence margin was preregistered.

The primary estimand of the underlying analysis is the **genus-level mean
difference**; every figure element resamples (or describes) genera, not
individual genomes.

Boundary: growth rates are codon-usage **predictions** of maximum growth
potential, not measured growth; candidate-carrier labels remain candidate-only
and do not establish a PHB/PHA degradation phenotype.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / ".mplconfig"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib

matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches
from matplotlib.collections import PolyCollection
from scipy import stats

from grodon_reanalysis_v2 import (
    CANDIDATE_CARRIER,
    CANDIDATE_NOT_DETECTED,
    GROWTH_RATE_UNIT,
    LEGACY_PHENOTYPE_LABELS,
    resolve_growth_rate,
)

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "runs" / "20260920_phaded_grodon_growth_01" / "results"
OUT_DIR = RUN / "figures" / "nature"
SOURCE_DIR = OUT_DIR / "source_data"
LABEL = "newdeg40k"
FIGURE_NAME = "growth_comparison_nature"

# ------------------------------------------------------------- text labels ---

#: Exact v2 vocabulary.  ``emit`` values are the only group labels this module
#: is allowed to put on an axis, in a legend or in a caption.
TEXT = {
    "primary_group_label": "Candidate carrier",
    "control_group_label": "Candidate not detected",
    "control_group_note": "not detected under the defined search",
    "secondary_tag": "Secondary",
    "primary_estimand_label": "Genus-level mean difference (primary estimand)",
    "resampling_unit_label": "resampling unit: genus",
    "per_genome_caption": "secondary (genomes are not independent); genus test in panel b",
    "compartment_caption": "secondary analysis; genus deltas, permutation test",
    "footnote": (
        "Growth rates are gRodon2 predictions of maximum growth potential, not "
        "measured growth; group labels are candidate-only (candidate gene carrier "
        "vs no candidate detected under the defined search) and do not establish a "
        "PHB/PHA degradation phenotype."
    ),
    "title": (
        "gRodon2 predicted maximum growth rate: candidate-carrier genomes vs "
        "same-genus genomes with no candidate detected under the defined search"
    ),
    "no_difference_wording": "no difference detected under the current design",
}

#: Label map by internal key, kept so callers never build a group label inline.
GROUP_LABELS = {
    CANDIDATE_CARRIER: TEXT["primary_group_label"],
    CANDIDATE_NOT_DETECTED: TEXT["control_group_label"],
}


def group_label(status: str) -> str:
    """Return the v2 display label for a candidate-detection status."""
    try:
        return GROUP_LABELS[status]
    except KeyError:
        raise ValueError(
            f"unknown candidate_detection_status {status!r}; expected one of "
            f"{sorted(GROUP_LABELS)}"
        ) from None


def secondary_label(text: str) -> str:
    """Prefix a secondary-analysis label with an explicit marker."""
    return f"{text} ({TEXT['secondary_tag'].lower()})"


# ---------------------------------------------------------------- style ----
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"

PALETTE = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "teal": "#42949E",
    "violet": "#9A4D8E",
    "red_strong": "#B64342",
    "neutral_light": "#D8D8D8",
    "neutral_mid": "#8F8F8F",
    "neutral_dark": "#4D4D4D",
    "neutral_black": "#272727",
    "bg_blue": "#EEF4FA",
    "bg_teal": "#EAF6F5",
    "bg_red": "#F9ECEA",
    "accent": "#2E8B57",
}

GROUP_COLOR = {"intracellular": PALETTE["blue_main"], "extracellular": PALETTE["teal"]}
GROUP_LABEL = {"intracellular": "Intracellular", "extracellular": "Extracellular"}

PANEL_LETTERS: list[tuple[str, object, object]] = []


def apply_style(font_size: float = 7.0) -> None:
    mpl.rcParams.update(
        {
            "font.size": font_size,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "savefig.dpi": 600,
        }
    )


def panel_label(ax, label: str) -> None:
    txt = ax.text(
        -0.16,
        1.20,
        label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
        ha="left",
        color=PALETTE["neutral_black"],
    )
    PANEL_LETTERS.append((label, ax, txt))


def caption_line(ax, text: str) -> None:
    """Small grey caption placed left-aligned between the axes and its title."""
    ax.text(
        0.0,
        1.008,
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=5.1,
        color=PALETTE["neutral_mid"],
    )


def title_artists(ax) -> list:
    """Every title artist of an axes, in reading order.

    ``set_title(..., loc="left")`` -- the convention used by every panel here --
    stores the string in the axes' *left-title* artist, not in ``ax.title``.
    Collecting only ``ax.title`` silently misses the panel titles.
    """
    artists = [ax.title]
    for attribute in ("_left_title", "_right_title"):
        artist = getattr(ax, attribute, None)
        if artist is not None and artist is not ax.title:
            artists.append(artist)
    return artists


def qc_panel_letters(fig) -> list[str]:
    """Return a list of panel letters that overlap their y-label or title."""
    fig.canvas.draw()
    problems = []
    for label, ax, txt in PANEL_LETTERS:
        bb = txt.get_window_extent()
        hits = []
        if ax.yaxis.label.get_text().strip() and bb.overlaps(ax.yaxis.label.get_window_extent()):
            hits.append("ylabel")
        for artist in title_artists(ax):
            if artist.get_text().strip() and bb.overlaps(artist.get_window_extent()):
                hits.append("title")
                break
        if hits:
            problems.append(f"{label}:{'+'.join(hits)}")
    return problems


def collect_figure_text(fig) -> list[str]:
    """Every piece of text one rendered figure emits (for label auditing).

    Covers titles, axis labels, tick labels, legends, annotations and figure
    text, so a test can assert on the *rendered* vocabulary rather than on the
    module source.
    """
    texts: list[str] = []
    for ax in fig.axes:
        candidates = title_artists(ax) + [ax.xaxis.label, ax.yaxis.label]
        candidates.extend(ax.texts)
        candidates.extend(ax.get_xticklabels())
        candidates.extend(ax.get_yticklabels())
        legend = ax.get_legend()
        if legend is not None:
            candidates.extend(legend.get_texts())
            if legend.get_title() is not None:
                candidates.append(legend.get_title())
        for artist in candidates:
            value = artist.get_text() if hasattr(artist, "get_text") else ""
            if value and value.strip():
                texts.append(value)
    texts.extend(t.get_text() for t in fig.texts if t.get_text().strip())
    if fig._suptitle is not None and fig._suptitle.get_text().strip():
        texts.append(fig._suptitle.get_text())
    return texts


def audit_figure_labels(fig) -> dict:
    """Check the rendered vocabulary for forbidden phenotype labels.

    ``degrader`` / ``non-degrader`` asserted a phenotype that a candidate
    classification does not establish; they must not come back through an axis,
    a legend, a title or a caption.
    """
    texts = collect_figure_text(fig)
    haystack = "\n".join(texts).lower()
    findings = [token for token in LEGACY_PHENOTYPE_LABELS if token in haystack]
    missing_vocabulary = [
        token
        for token in (TEXT["primary_group_label"], TEXT["control_group_label"])
        if token.lower() not in haystack
    ]
    report = {
        "schema_version": 1,
        "texts_measured": len(texts),
        "forbidden_labels": list(LEGACY_PHENOTYPE_LABELS),
        "forbidden_labels_found": findings,
        "required_vocabulary_missing": missing_vocabulary,
        "verdict": "PASS" if not findings else "FIX BEFORE DELIVERY",
    }
    return report


def assert_figure_labels(fig) -> None:
    """Raise ``ValueError`` when a rendered figure emits a forbidden label."""
    report = audit_figure_labels(fig)
    if report["forbidden_labels_found"]:
        raise ValueError(
            "forbidden legacy phenotype label in figure text: "
            f"{report['forbidden_labels_found']}"
        )


def save_figure(fig, name: str, out_dir: Path | None = None) -> None:
    """Export the vector pair, a 600-dpi TIFF submission raster and a PNG preview."""
    out_dir = OUT_DIR if out_dir is None else Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{name}.svg", bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.tiff", dpi=600, bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def audit_render_collisions(fig, name: str, out_dir: Path | None = None) -> dict:
    """Backend-native rendered collision audit.

    Substitute for the skill's ``audit_figure_collisions.py``, which needs
    PyMuPDF (unavailable here). It covers the same two FAIL classes from the
    same final rendered figure, measured with the drawing backend's own
    geometry: ``text-text`` (two text boxes materially overlap) and
    ``text-stroke`` (a stroked line/curve point falls inside a text box).
    """
    from matplotlib.transforms import Bbox

    out_dir = OUT_DIR if out_dir is None else Path(out_dir)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    texts = []
    for ax in fig.axes:
        for artist in list(ax.texts) + [ax.title, ax.xaxis.label, ax.yaxis.label]:
            if artist.get_text().strip() and artist.get_visible():
                texts.append(artist)
        texts.extend(t for t in ax.get_xticklabels() if t.get_text().strip())
        texts.extend(t for t in ax.get_yticklabels() if t.get_text().strip())
    texts.extend(t for t in fig.texts if t.get_text().strip())

    boxes = [(t, t.get_window_extent(renderer)) for t in texts]

    def where(artist) -> str:
        owner = getattr(artist, "axes", None)
        if owner is not None and owner in fig.axes:
            return f"axes{fig.axes.index(owner)}"
        return "figure"

    text_text = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            first, second = boxes[i][1], boxes[j][1]
            overlap = Bbox.intersection(first, second)
            if overlap is not None and overlap.width > 1.0 and overlap.height > 1.0:
                text_text.append(
                    {
                        "a": boxes[i][0].get_text()[:40],
                        "b": boxes[j][0].get_text()[:40],
                        "where_a": where(boxes[i][0]),
                        "where_b": where(boxes[j][0]),
                        "overlap_pt": [round(overlap.width, 2), round(overlap.height, 2)],
                    }
                )

    text_stroke = []
    for ax in fig.axes:
        for line in ax.lines:
            if not line.get_visible():
                continue
            points = line.get_xydata()
            if len(points) == 0:
                continue
            try:
                # Use the artist's OWN transform: axvline/axhline live in
                # blended (xaxis/yaxis) coordinates, not ax.transData.
                pixels = line.get_transform().transform(np.asarray(points, dtype=float))
            except Exception:
                continue
            for label, box in boxes:
                x0, y0, x1, y1 = box.extents
                hit = 0
                # Sample ALONG each segment: a straight line can cross a text box
                # with both endpoints outside it (a full-height axvline does), so
                # vertex-only testing silently misses that class.
                for k in range(len(pixels) - 1):
                    segment = np.linspace(pixels[k], pixels[k + 1], 64)
                    inside = (
                        (segment[:, 0] >= x0)
                        & (segment[:, 0] <= x1)
                        & (segment[:, 1] >= y0)
                        & (segment[:, 1] <= y1)
                    )
                    hit += int(inside.sum())
                if hit:
                    text_stroke.append({"text": label.get_text()[:40], "line_points_inside": hit})
                    break

    # text-fill-edge WARN class: text partly sitting on bars / filled violin bodies.
    fill_artists = []
    for ax in fig.axes:
        for patch in ax.patches:
            try:
                if patch.get_width() * patch.get_height() > 0:
                    fill_artists.append((ax, patch))
            except AttributeError:
                continue
        for coll in ax.collections:
            if isinstance(coll, PolyCollection) and len(coll.get_paths()) > 0:
                fill_artists.append((ax, coll))

    text_fill_edge = []
    text_fill_contained = []
    for ax, artist in fill_artists:
        try:
            fill_box = artist.get_window_extent(renderer)
        except Exception:
            continue
        # Only the part of a fill that is actually visible inside its axes can
        # collide with text (PDF-geometry equivalent: renderers clip at the axes).
        axes_box = ax.get_window_extent(renderer)
        clipped = Bbox.intersection(fill_box, axes_box)
        if clipped is None:
            continue
        fill_box = clipped
        if fill_box.width < 2 or fill_box.height < 2:
            continue
        for label, box in boxes:
            if box.x1 < fill_box.x0 or box.x0 > fill_box.x1 or box.y1 < fill_box.y0 or box.y0 > fill_box.y1:
                continue
            overlap_x = min(box.x1, fill_box.x1) - max(box.x0, fill_box.x0)
            overlap_y = min(box.y1, fill_box.y1) - max(box.y0, fill_box.y0)
            if overlap_x <= 1.0 or overlap_y <= 1.0:
                continue
            entry = {
                "text": label.get_text()[:40],
                "where": where(label),
                "fill": type(artist).__name__,
                "overlap_pt": [round(overlap_x, 2), round(overlap_y, 2)],
            }
            contained = (
                box.x0 >= fill_box.x0
                and box.x1 <= fill_box.x1
                and box.y0 >= fill_box.y0
                and box.y1 <= fill_box.y1
            )
            (text_fill_contained if contained else text_fill_edge).append(entry)

    # crowding check: text sitting within CLEARANCE_PT of ANOTHER panel's plot area
    # (not a collision class of the contract, but the gutter case that hides behind it).
    clearance_pt = 2.0
    text_axes_crowding = []
    axes_rects = [(ax, ax.get_window_extent(renderer)) for ax in fig.axes]
    for label, box in boxes:
        owner = getattr(label, "axes", None)
        for ax, rect in axes_rects:
            if ax is owner:
                continue
            if box.x1 < rect.x0 or box.x0 > rect.x1 or box.y1 < rect.y0 or box.y0 > rect.y1:
                gap_x = max(rect.x0 - box.x1, box.x0 - rect.x1, 0.0)
                gap_y = max(rect.y0 - box.y1, box.y0 - rect.y1, 0.0)
                if max(gap_x, gap_y) <= clearance_pt:
                    text_axes_crowding.append(
                        {
                            "text": label.get_text()[:40],
                            "where": where(label),
                            "other_axes": f"axes{fig.axes.index(ax)}",
                            "gap_pt": [round(gap_x, 2), round(gap_y, 2)],
                        }
                    )

    # clipping check: a patch whose geometry leaves its axes is silently cut
    # (e.g. a rounded box drawn flush to the axes edge loses its border).
    clipped_artists = []
    for ax in fig.axes:
        axes_box = ax.get_window_extent(renderer)
        for patch in ax.patches:
            try:
                extent = patch.get_window_extent(renderer)
            except Exception:
                continue
            if extent.width <= 0 or extent.height <= 0:
                continue
            overflow_x = max(axes_box.x0 - extent.x0, extent.x1 - axes_box.x1, 0.0)
            overflow_y = max(axes_box.y0 - extent.y0, extent.y1 - axes_box.y1, 0.0)
            if overflow_x > 0.2 or overflow_y > 0.2:
                clipped_artists.append(
                    {
                        "artist": type(patch).__name__,
                        "where": f"axes{fig.axes.index(ax)}",
                        "overflow_pt": [round(overflow_x, 2), round(overflow_y, 2)],
                    }
                )

    report = {
        "schema_version": 1,
        "backend": "matplotlib",
        "method": "backend_native_geometry",
        "substitute_for": "nature-figure audit_figure_collisions.py (PyMuPDF unavailable)",
        "texts_measured": len(boxes),
        "text_text": text_text,
        "text_stroke": text_stroke,
        "text_fill_edge": text_fill_edge,
        "text_fill_contained": text_fill_contained,
        "text_axes_crowding": text_axes_crowding,
        "clipped_artists": clipped_artists,
        "verdict": "FIX BEFORE DELIVERY"
        if (text_text or text_stroke or text_fill_edge or text_axes_crowding or clipped_artists)
        else "PASS",
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{name}.collision-audit.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"[collision] verdict={report['verdict']} texts={len(boxes)} "
        f"text-text={len(text_text)} text-stroke={len(text_stroke)} "
        f"text-fill-edge={len(text_fill_edge)} contained={len(text_fill_contained)} "
        f"crowding={len(text_axes_crowding)} clipped={len(clipped_artists)}"
    )
    return report


def run_alignment_gate(fig, name: str, out_dir: Path | None = None) -> None:
    """Mandatory render-time multi-panel alignment gate (nature-figure contract)."""
    out_dir = OUT_DIR if out_dir is None else Path(out_dir)
    scripts = ROOT / ".superpowers" / "nature-skills" / "skills" / "nature-figure" / "scripts"
    if not (scripts / "audit_panel_alignment.py").exists():
        print("[alignment] nature-figure skill scripts not found; gate not run (documented exception)")
        return
    sys.path.insert(0, str(scripts))
    from audit_panel_alignment import require_matplotlib_panel_alignment

    report = require_matplotlib_panel_alignment(
        fig,
        json_out=str(out_dir / f"{name}.alignment.json"),
        overlay_svg=str(out_dir / f"{name}.alignment.svg"),
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        require_panel_labels=True,
        strict=True,
    )
    print(f"[alignment] verdict={report.get('verdict')} panels={len((report.get('layout') or {}).get('panels', []))}")


# ------------------------------------------------------------- helpers ----
def cliff_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta for a vs b (positive => a tends to be larger)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    gt = sum(float(np.sum(ai > b)) for ai in a)
    lt = sum(float(np.sum(ai < b)) for ai in a)
    return (gt - lt) / (len(a) * len(b))


def genus_bootstrap_ci(
    a: np.ndarray, b: np.ndarray, rng: np.random.Generator, n_boot: int = 2000
) -> tuple[float, float]:
    """Bootstrap 95% CI of mean(b) - mean(a) for one genus."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) == 0 or len(b) == 0:
        return (float("nan"), float("nan"))
    est = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        est[i] = (
            b[rng.integers(0, len(b), size=len(b))].mean()
            - a[rng.integers(0, len(a), size=len(a))].mean()
        )
    lo, hi = np.percentile(est, [2.5, 97.5])
    return float(lo), float(hi)


def metric(tests: pd.DataFrame, name: str) -> float:
    value = tests.loc[tests["metric"] == name, "value"]
    if value.empty:
        raise KeyError(f"missing metric: {name}")
    return float(value.iloc[0])


def bootstrap_mean_ci(values: np.ndarray, rng: np.random.Generator, n_boot: int = 1200) -> tuple[float, float]:
    """Bootstrap 95% CI of the mean of ``values`` (resampling the values)."""
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return (float("nan"), float("nan"))
    est = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        est[i] = values[rng.integers(0, len(values), size=len(values))].mean()
    lo, hi = np.percentile(est, [2.5, 97.5])
    return float(lo), float(hi)


def verify_growth_rate_unit(balanced: pd.DataFrame) -> dict:
    """Verify the resolved growth-rate column against the declared formula.

    The balanced table carries ``growth_rate_per_h``; when the source
    doubling-time column is available as well, exactly one conversion is
    expected (``ln2 / d``).  This makes a wrong or double conversion a hard
    failure instead of a silently mis-scaled axis.
    """
    columns = set(balanced.columns)
    rate_column = "growth_rate_per_h" if "growth_rate_per_h" in columns else None
    if rate_column is None:
        raise ValueError(
            "the balanced table has no growth_rate_per_h column; the growth-rate "
            "unit cannot be verified"
        )
    report = {
        "rate_column": rate_column,
        "unit": GROWTH_RATE_UNIT,
        "formula": "growth_rate_per_h = ln2 / doubling_time_h",
        "checked_rows": 0,
        "max_relative_error": 0.0,
        "conversions_applied": "one (doubling time -> rate); the rate column is never converted again",
    }
    if "doubling_time_h" not in columns:
        report["cross_check"] = "not available: no doubling_time_h column"
        return report
    frame = balanced[["doubling_time_h", rate_column]].apply(
        pd.to_numeric, errors="coerce"
    ).dropna()
    if frame.empty:
        report["cross_check"] = "not available: no numeric (doubling_time_h, rate) pair"
        return report
    derived = np.log(2.0) / frame["doubling_time_h"].to_numpy(float)
    recorded = frame[rate_column].to_numpy(float)
    relative = np.abs(recorded - derived) / derived
    report["checked_rows"] = int(len(frame))
    report["max_relative_error"] = float(np.max(relative))
    if report["max_relative_error"] > 0.05:
        raise ValueError(
            "growth_rate_per_h does not match ln2/doubling_time_h "
            f"(max relative error {report['max_relative_error']:.3g}); the unit or "
            "the conversion count is wrong"
        )
    report["cross_check"] = "consistent with ln2/doubling_time_h"
    return report


class FigureInputs:
    """Everything the six panels need; loaded from disk or supplied by a test.

    A plain class rather than a ``dataclass``: this module is loaded by tests
    through ``importlib.util.spec_from_file_location`` without being registered
    in ``sys.modules``, and ``dataclass`` resolves annotations through
    ``sys.modules``.
    """

    __slots__ = (
        "balanced",
        "summary",
        "tests",
        "sf_tests",
        "sf_effects",
        "grp_tests",
        "grp_effects",
        "manifest_stats",
        "dedup_stats",
    )

    def __init__(
        self,
        *,
        balanced: pd.DataFrame,
        summary: pd.DataFrame,
        tests: pd.DataFrame,
        sf_tests: pd.DataFrame,
        sf_effects: pd.DataFrame,
        grp_tests: pd.DataFrame,
        grp_effects: pd.DataFrame,
        manifest_stats: dict,
        dedup_stats: dict,
    ) -> None:
        self.balanced = balanced
        self.summary = summary
        self.tests = tests
        self.sf_tests = sf_tests
        self.sf_effects = sf_effects
        self.grp_tests = grp_tests
        self.grp_effects = grp_effects
        self.manifest_stats = manifest_stats
        self.dedup_stats = dedup_stats


def load_inputs() -> FigureInputs:
    balanced = pd.read_csv(RUN / f"grodon_growth_balanced_by_genus_{LABEL}.tsv", sep="\t")
    balanced["growth_rate_per_h"] = pd.to_numeric(balanced["growth_rate_per_h"], errors="coerce")
    balanced = balanced[balanced["growth_rate_per_h"].notna()].copy()
    return FigureInputs(
        balanced=balanced,
        summary=pd.read_csv(RUN / f"grodon_growth_balanced_by_genus_summary_{LABEL}.tsv", sep="\t"),
        tests=pd.read_csv(RUN / f"grodon_growth_statistical_tests_{LABEL}.tsv", sep="\t"),
        sf_tests=pd.read_csv(RUN / f"grodon_growth_superfamily_statistical_tests_{LABEL}.tsv", sep="\t"),
        sf_effects=pd.read_csv(RUN / f"grodon_growth_superfamily_effects_{LABEL}.tsv", sep="\t"),
        grp_tests=pd.read_csv(RUN / f"grodon_growth_group_statistical_tests_{LABEL}.tsv", sep="\t"),
        grp_effects=pd.read_csv(RUN / f"grodon_growth_group_effects_{LABEL}.tsv", sep="\t"),
        manifest_stats=pd.read_json(RUN / "manifest_stats.json", typ="series").to_dict(),
        dedup_stats=pd.read_json(RUN / "predictions_dedup_stats.json", typ="series").to_dict(),
    )


# ---------------------------------------------------------------- main ----
def build_figure(
    data: FigureInputs,
    *,
    rng: np.random.Generator | None = None,
    out_dir: Path | None = None,
    source_dir: Path | None = None,
    write_audits: bool = True,
) -> object:
    """Draw the six-panel figure and (optionally) write audits and source data.

    Split out of :func:`main` so the rendered vocabulary can be audited from a
    synthetic table without reading ``runs/`` or writing a real figure.
    """
    apply_style()
    out_dir = OUT_DIR if out_dir is None else Path(out_dir)
    source_dir = SOURCE_DIR if source_dir is None else Path(source_dir)
    rng = np.random.default_rng(np.random.PCG64(20260920)) if rng is None else rng

    balanced = data.balanced
    summary = data.summary
    tests = data.tests
    sf_tests = data.sf_tests
    sf_effects = data.sf_effects
    grp_tests = data.grp_tests
    grp_effects = data.grp_effects
    manifest_stats = data.manifest_stats
    dedup_stats = data.dedup_stats

    unit_report = verify_growth_rate_unit(balanced)
    print(f"[units] {unit_report}")

    pos = balanced[balanced["candidate_detection_status"] == CANDIDATE_CARRIER]
    neg = balanced[balanced["candidate_detection_status"] == CANDIDATE_NOT_DETECTED]
    g_pos = pos["growth_rate_per_h"].to_numpy(float)
    g_neg = neg["growth_rate_per_h"].to_numpy(float)

    n_pos = int(metric(tests, "n_balanced_candidate_carrier"))
    n_neg = int(metric(tests, "n_balanced_candidate_not_detected"))
    n_genera = int(metric(tests, "n_balanced_genera"))
    mean_delta = metric(tests, "genus_mean_delta_growth_rate_per_h")
    median_delta = metric(tests, "genus_median_delta_growth_rate_per_h")
    p_wilcox = metric(tests, "wilcoxon_signed_rank_two_sided_p")
    p_perm = metric(tests, "stratified_permutation_unweighted_p")
    p_sign = metric(tests, "exact_sign_test_two_sided_p")
    effect_r = metric(tests, "wilcoxon_effect_r")

    fig = plt.figure(figsize=(7.2, 5.6))
    # Vertical/horizontal gutters are wider than in the 2026-09-20 revision of
    # this figure: the v2 group labels occupy two wrapped lines, and the extra
    # row height keeps panel a's annotation clear of its own tick labels.
    gs = fig.add_gridspec(
        2, 3, height_ratios=[1.0, 0.92], width_ratios=[1.0, 1.0, 1.0], hspace=0.88, wspace=1.50
    )

    # ---------------------------------------------------- panel a --------
    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "a")
    bp = ax.boxplot(
        [g_neg, g_pos],
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "lw": 1.0},
        boxprops={"linewidth": 0.7},
        whiskerprops={"linewidth": 0.7},
        capprops={"linewidth": 0.7},
    )
    bp["boxes"][0].set_facecolor(PALETTE["bg_red"])
    bp["boxes"][0].set_edgecolor(PALETTE["red_strong"])
    bp["boxes"][1].set_facecolor(PALETTE["bg_blue"])
    bp["boxes"][1].set_edgecolor(PALETTE["blue_main"])
    ax.scatter(1 + (rng.random(len(g_neg)) - 0.5) * 0.11, g_neg, s=4, color=PALETTE["red_strong"], alpha=0.28, edgecolors="none", zorder=3)
    ax.scatter(2 + (rng.random(len(g_pos)) - 0.5) * 0.11, g_pos, s=4, color=PALETTE["blue_main"], alpha=0.28, edgecolors="none", zorder=3)
    p_mw = stats.mannwhitneyu(g_pos, g_neg, alternative="two-sided").pvalue
    cd = cliff_delta(g_pos, g_neg)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(
        [
            f"{TEXT['control_group_label']}\nn = {n_neg:,}",
            f"{TEXT['primary_group_label']}\nn = {n_pos:,}",
        ],
        fontsize=6.2,
    )
    ax.set_ylabel(f"Predicted max growth rate ({GROWTH_RATE_UNIT})")
    ax.set_title(secondary_label("Genome-level distributions"), loc="left", fontsize=7.6, pad=10)
    ax.text(
        0.03,
        0.72,
        f"Mann-Whitney P = {p_mw:.2g}\nCliff's δ = {cd:+.3f}\n"
        f"medians {np.median(g_neg):.3f} vs {np.median(g_pos):.3f} ({GROWTH_RATE_UNIT})",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.2,
    )
    caption_line(ax, TEXT["per_genome_caption"])

    # ---------------------------------------------------- panel b --------
    ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "b")
    deltas = summary["delta_carrier_minus_control"].to_numpy(float)
    lo, hi = np.percentile(deltas, [1, 99])
    bins = np.linspace(lo, hi, 46)
    ax.hist(np.clip(deltas, lo, hi), bins=bins, color="#9ECAE1", edgecolor="white", linewidth=0.3)
    ax.axvline(0, color=PALETTE["neutral_mid"], lw=0.8, ls="--", zorder=1, ymax=0.62)
    ax.axvline(median_delta, color=PALETTE["blue_main"], lw=1.3, ymax=0.62)
    ax.axvline(mean_delta, color=PALETTE["red_strong"], lw=1.3, ls="-", ymax=0.62)
    ax.set_xlabel(f"Genus-level delta ({GROWTH_RATE_UNIT})")
    ax.set_ylabel("Genera")
    ax.set_ylim(0, ax.get_ylim()[1] * 1.62)
    ax.set_title(TEXT["primary_estimand_label"], loc="left", fontsize=7.6, pad=10)
    ax.text(
        0.98,
        0.96,
        f"n = {n_genera:,} genera\nmean = {mean_delta:+.4f}\nmedian = {median_delta:+.4f}\n"
        f"signed-rank P = {p_wilcox:.3g}\npermutation P = {p_perm:.3g}\nsign test P = {p_sign:.3g}\n"
        f"effect r = {effect_r:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.9,
    )
    caption_line(ax, "lines: mean (red), median (blue); x: 1st-99th pct; " + TEXT["resampling_unit_label"])

    # ---------------------------------------------------- panel c --------
    ax = fig.add_subplot(gs[0, 2])
    panel_label(ax, "c")
    sizes = 10 + 1.6 * np.sqrt(summary["n_candidate_carrier"].to_numpy(float))
    sc = ax.scatter(
        summary["mean_growth_candidate_not_detected"],
        summary["mean_growth_candidate_carrier"],
        c=deltas,
        cmap="RdBu_r",
        s=sizes,
        alpha=0.65,
        edgecolor="white",
        linewidth=0.2,
        vmin=-np.percentile(np.abs(deltas), 98),
        vmax=np.percentile(np.abs(deltas), 98),
    )
    lim = max(summary["mean_growth_candidate_not_detected"].max(), summary["mean_growth_candidate_carrier"].max()) * 1.22
    ax.plot([0, lim], [0, lim], color=PALETTE["neutral_mid"], lw=0.8, ls="--", zorder=1)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel(f"{TEXT['control_group_label']} mean ({GROWTH_RATE_UNIT})")
    ax.set_ylabel(f"{TEXT['primary_group_label']} mean ({GROWTH_RATE_UNIT})")
    ax.set_title("Paired genus means", loc="left", fontsize=7.6, pad=10)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label(f"delta ({GROWTH_RATE_UNIT})", fontsize=6.2)
    cbar.ax.tick_params(labelsize=5.6)
    above = int((deltas > 0).sum())
    caption_line(
        ax,
        f"{above:,}/{n_genera:,} genera above 1:1; size = pairs; {TEXT['resampling_unit_label']}",
    )

    # ---------------------------------------------------- panel d --------
    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "d")
    intra = grp_effects.loc[grp_effects["group"] == "intracellular", "delta_carrier_minus_control"].to_numpy(float)
    extra = grp_effects.loc[grp_effects["group"] == "extracellular", "delta_carrier_minus_control"].to_numpy(float)
    parts = ax.violinplot([intra, extra], widths=0.72, showextrema=False)
    for body, key in zip(parts["bodies"], ("intracellular", "extracellular")):
        body.set_facecolor(PALETTE["bg_blue"] if key == "intracellular" else PALETTE["bg_teal"])
        body.set_edgecolor(GROUP_COLOR[key])
        body.set_alpha(0.95)
        body.set_linewidth(0.7)
    for i, (vals, key) in enumerate(((intra, "intracellular"), (extra, "extracellular")), start=1):
        ax.scatter(i + (rng.random(len(vals)) - 0.5) * 0.09, vals, s=3, color=GROUP_COLOR[key], alpha=0.22, edgecolors="none", zorder=3)
        ax.plot([i - 0.28, i + 0.28], [np.median(vals)] * 2, color=PALETTE["neutral_black"], lw=1.1, zorder=4)
        ax.scatter([i], [vals.mean()], marker="D", s=18, color=PALETTE["red_strong"], zorder=5)
    row = grp_tests.loc[grp_tests["group"] == "intracellular_vs_extracellular"].iloc[0]
    ax.set_xticks([1, 2])
    ax.set_xticklabels(
        [
            f"{GROUP_LABEL['intracellular']}\nn = {len(intra)} genera",
            f"{GROUP_LABEL['extracellular']}\nn = {len(extra)} genera",
        ],
        fontsize=6.2,
    )
    ax.set_ylabel(f"Genus-level delta ({GROWTH_RATE_UNIT})")
    ax.set_title(secondary_label("Intracellular vs extracellular"), loc="left", fontsize=7.6, pad=10)
    all_d = np.concatenate([intra, extra])
    d_lo, d_hi = (float(v) for v in np.percentile(all_d, [2, 98]))
    ax.set_ylim(d_lo, d_hi)
    caption_line(
        ax,
        f"{TEXT['compartment_caption']}; Δ = {float(row['mean_difference_intra_minus_extra']):+.4f} "
        f"({GROWTH_RATE_UNIT}); Mann-Whitney P = {float(row['mannwhitney_u_p']):.2g}",
    )

    # ---------------------------------------------------- panel e --------
    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "e")
    sf_order = sf_tests.sort_values("mean_delta", ascending=True)["subtype"].tolist()
    ys = np.arange(len(sf_order))
    for y, subtype in zip(ys, sf_order):
        vals = sf_effects.loc[sf_effects["subtype"] == subtype, "delta_carrier_minus_control"].to_numpy(float)
        lo_ci, hi_ci = bootstrap_mean_ci(vals, rng)
        mean_v = float(np.mean(vals))
        q_row = sf_tests.loc[sf_tests["subtype"] == subtype].iloc[0]
        ax.plot([lo_ci, hi_ci], [y, y], color=PALETTE["neutral_mid"], lw=1.5, alpha=0.8, solid_capstyle="butt", zorder=2)
        ax.scatter([mean_v], [y], s=20, color=PALETTE["blue_main"] if mean_v > 0 else PALETTE["red_strong"], zorder=3, edgecolors="white", linewidth=0.3)
        ax.text(hi_ci + 0.02, y, f"  q = {float(q_row['wilcoxon_q_bh']):.2f}", va="center", ha="left", fontsize=5.3, color=PALETTE["neutral_mid"])
    ax.axvline(0, color=PALETTE["neutral_mid"], lw=0.8, ls="--", zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels(
        [
            s.replace("extracellular dPHASCL type 1", "Extra. dPHASCL type 1")
            .replace("extracellular dPHASCL type 2", "Extra. dPHASCL type 2")
            .replace("extracellular native-SCL/PhaZ7-like", "Extra. PhaZ7-like")
            .replace("extracellular dPHAMCL", "Extra. dPHAMCL")
            .replace("intracellular nPHASCL without lipase box", "Intra. nPHASCL\n(no lipase box)")
            .replace("intracellular nPHAMCL", "Intra. nPHAMCL")
            for s in sf_order
        ],
        fontsize=5.4,
        linespacing=1.2,
    )
    ax.set_xlabel(f"Mean genus delta ({GROWTH_RATE_UNIT})")
    ax.set_title("Superfamily-stratified", loc="left", fontsize=7.6, pad=10)
    x_lo = min(-0.05, float(sf_tests["mean_delta"].min()) - 0.10)
    x_hi = float(sf_tests["mean_delta"].max()) + 0.30
    ax.set_xlim(x_lo, x_hi)
    caption_line(ax, "error bars: bootstrap 95% CI over genera; q = BH")

    # ---------------------------------------------------- panel f --------
    ax = fig.add_subplot(gs[1, 2])
    panel_label(ax, "f")
    ax.axis("off")
    boxes = [
        ("1  Annotations", "40,690"),
        ("2  Genomes", "30,623"),
        ("3  Manifest", f"{manifest_stats['manifest_carrier']:,}+{manifest_stats['manifest_control']:,}"),
        ("4  Analysis", f"{n_pos + n_neg:,}"),
    ]
    positions = [(0.022, 0.785), (0.520, 0.785), (0.022, 0.585), (0.520, 0.585)]
    for (title, value), (x, y) in zip(boxes, positions):
        ax.add_patch(
            patches.FancyBboxPatch(
                (x, y),
                0.456,
                0.175,
                boxstyle="round,pad=0.008,rounding_size=0.014",
                facecolor="#F5F5F5",
                edgecolor="#BDBDBD",
                linewidth=0.7,
                transform=ax.transAxes,
            )
        )
        ax.text(x + 0.018, y + 0.118, title, ha="left", va="center", fontsize=5.3)
        ax.text(x + 0.018, y + 0.048, value, ha="left", va="center", fontsize=6.4, fontweight="bold")
    ax.text(
        0.01,
        0.54,
        f"status = ok: {int(dedup_stats['ok']):,};  failed: {int(dedup_stats['failed'])} "
        f"({100 * dedup_stats['failed'] / dedup_stats['unique_genomes']:.2f}%)\n"
        f"no same-genus control: {manifest_stats['carrier_genomes_without_a_same_genus_control']:,} "
        f"({100 * manifest_stats['carrier_genomes_without_a_same_genus_control'] / manifest_stats['carrier_genomes_input']:.1f}%)\n"
        f"manifest genera: {manifest_stats['manifest_genera']:,};  analysed: {n_genera:,}\n"
        f"{TEXT['primary_estimand_label']}\n"
        f"Carriers: small positive shift (permutation P = {p_perm:.3g}; r = {effect_r:.2f})\n"
        f"Intra vs extra: {TEXT['no_difference_wording']} (P = {float(row['p_value']):.2f})\n"
        "Predictions; candidate-only labels (no phenotype claim).",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.4,
        color=PALETTE["neutral_black"],
        linespacing=1.5,
    )

    fig.suptitle(TEXT["title"], fontsize=8.6, y=0.985)
    fig.text(
        0.5,
        0.008,
        f"GTDB R232; genus-balanced design ({n_genera:,} genera); panels b and d axes clipped to the 2nd-98th percentile. "
        + TEXT["footnote"],
        ha="center",
        va="bottom",
        fontsize=5.6,
        color=PALETTE["neutral_mid"],
    )
    fig.subplots_adjust(left=0.075, right=0.975, top=0.885, bottom=0.105)

    label_report = audit_figure_labels(fig)
    print(f"[labels] {label_report}")
    assert_figure_labels(fig)

    problems = qc_panel_letters(fig)
    print(f"[QC] panel-letter overlaps: {problems if problems else 'none'}")
    if write_audits:
        run_alignment_gate(fig, FIGURE_NAME, out_dir)
        audit_render_collisions(fig, FIGURE_NAME, out_dir)
        save_figure(fig, FIGURE_NAME, out_dir)

        # --------------------------------------------------- source data -
        source_dir.mkdir(parents=True, exist_ok=True)
        balanced[
            [
                "genome_id",
                "candidate_detection_status",
                "group",
                "major_subtype",
                "phylum",
                "genus",
                "species",
                "growth_rate_per_h",
                "doubling_time_h",
            ]
        ].to_csv(source_dir / "source_per_genome.tsv", sep="\t", index=False)
        summary.to_csv(source_dir / "source_genus_deltas.tsv", sep="\t", index=False)
        tests.to_csv(source_dir / "source_main_tests.tsv", sep="\t", index=False)
        sf_tests.to_csv(source_dir / "source_superfamily_tests.tsv", sep="\t", index=False)
        grp_tests.to_csv(source_dir / "source_group_tests.tsv", sep="\t", index=False)
        pd.DataFrame(
            [
                {"stage": "candidate annotations", "value": 40690},
                {"stage": "genomes mapped", "value": manifest_stats["carrier_genomes_input"]},
                {"stage": "same-genus manifest", "value": manifest_stats["manifest_rows"]},
                {"stage": "predicted ok", "value": int(dedup_stats["ok"])},
                {"stage": "balanced analysis", "value": n_pos + n_neg},
            ]
        ).to_csv(source_dir / "source_ledger.tsv", sep="\t", index=False)

        print(f"wrote figure to {out_dir}")
        print(f"wrote source data to {source_dir}")
    return fig


def main() -> None:
    build_figure(load_inputs())


if __name__ == "__main__":
    main()
