#!/usr/bin/env python3
"""Generate SVG figures and summary tables for the LERW experiments."""

from __future__ import annotations

import argparse
import csv
import html
import math
import random
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, variance

from analyse_winding_results import GroupStats, grouped_stats, linear_fit, read_rows
from lerw_random_environment import (
    RandomEnvironment,
    loop_erase,
    random_walk_to_boundary,
    winding,
)


@dataclass(frozen=True)
class Series:
    label: str
    points: list[tuple[float, float]]
    color: str


COLORS = {
    "symmetric": "#1f2937",
    "gamma": "#dc2626",
    "gamma4": "#2563eb",
    "k20": "#0f766e",
    "k10": "#16a34a",
    "k5": "#84cc16",
    "k2": "#f59e0b",
    "k1": "#dc2626",
    "k0.5": "#7c2d12",
    "gamma4_k20": "#0f766e",
    "gamma4_k10": "#0891b2",
    "gamma4_k5": "#2563eb",
    "gamma4_k2": "#7c3aed",
    "gamma4_k1": "#db2777",
    "gamma4_k0.5": "#581c87",
}


def fmt_k(k: str | float | None) -> str:
    if k in {None, "", "None"}:
        return "None"
    value = float(k)
    return f"{value:g}"


def model_label(model: str, k: str | float | None) -> str:
    if model == "symmetric":
        return "symmetric"
    return f"{model}, k={fmt_k(k)}"


def key_color(model: str, k: str | float | None) -> str:
    if model == "symmetric":
        return COLORS["symmetric"]
    k_label = f"k{fmt_k(k)}"
    if model == "gamma4":
        return COLORS.get(f"gamma4_{k_label}", COLORS["gamma4"])
    return COLORS.get(k_label, COLORS["gamma"])


def stat_lookup(stats: list[GroupStats]) -> dict[tuple[str, str, int], GroupStats]:
    return {(item.model, item.k, item.L): item for item in stats}


def nice_ticks(v_min: float, v_max: float, count: int = 6) -> list[float]:
    if not math.isfinite(v_min) or not math.isfinite(v_max):
        return []
    if v_min == v_max:
        return [v_min]

    span = v_max - v_min
    rough = span / max(count - 1, 1)
    magnitude = 10 ** math.floor(math.log10(abs(rough)))
    fraction = rough / magnitude
    if fraction <= 1:
        nice_fraction = 1
    elif fraction <= 2:
        nice_fraction = 2
    elif fraction <= 5:
        nice_fraction = 5
    else:
        nice_fraction = 10

    step = nice_fraction * magnitude
    start = math.floor(v_min / step) * step
    end = math.ceil(v_max / step) * step
    ticks = []
    value = start
    while value <= end + step * 0.5:
        ticks.append(0.0 if abs(value) < 1e-12 else value)
        value += step
    return ticks


def format_tick(value: float) -> str:
    if abs(value) >= 1000:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    return f"{value:.2f}".rstrip("0").rstrip(".")


def svg_text(
    x: float,
    y: float,
    text: str,
    size: int = 13,
    anchor: str = "middle",
    weight: str = "400",
    fill: str = "#111827",
) -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" '
        f'fill="{fill}">{html.escape(text)}</text>'
    )


def svg_doc(width: int, height: int, body: list[str]) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        '<rect width="100%" height="100%" fill="#ffffff"/>\n'
        + "\n".join(body)
        + "\n</svg>\n"
    )


def write_svg(path: Path, width: int, height: int, body: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(svg_doc(width, height, body), encoding="utf-8")


def line_plot(
    path: Path,
    title: str,
    xlabel: str,
    ylabel: str,
    series_list: list[Series],
    x_ticks: list[tuple[float, str]] | None = None,
    y_min: float | None = None,
    y_max: float | None = None,
    width: int = 1120,
    height: int = 1120,
) -> None:
    if not series_list:
        return

    left, right, top, bottom = 88, 300, 58, 82
    plot_w = width - left - right
    plot_h = height - top - bottom

    xs = [x for series in series_list for x, _ in series.points]
    ys = [y for series in series_list for _, y in series.points]
    x_min, x_max = min(xs), max(xs)
    if x_min == x_max:
        x_min -= 1
        x_max += 1
    if y_min is None:
        y_min = min(ys)
    if y_max is None:
        y_max = max(ys)
    y_pad = (y_max - y_min) * 0.08 if y_max != y_min else 1.0
    y_min -= y_pad
    y_max += y_pad

    def px(x: float) -> float:
        return left + (x - x_min) / (x_max - x_min) * plot_w

    def py(y: float) -> float:
        return top + (y_max - y) / (y_max - y_min) * plot_h

    body: list[str] = []
    body.append(svg_text(width / 2, 28, title, size=19, weight="700"))
    body.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#f9fafb" stroke="#d1d5db"/>')

    y_ticks = nice_ticks(y_min, y_max)
    for tick in y_ticks:
        y = py(tick)
        body.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}" stroke="#e5e7eb"/>')
        body.append(svg_text(left - 10, y + 4, format_tick(tick), size=11, anchor="end", fill="#4b5563"))

    if x_ticks is None:
        x_ticks = [(tick, format_tick(tick)) for tick in nice_ticks(x_min, x_max)]
    for tick, label in x_ticks:
        x = px(tick)
        body.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_h}" stroke="#eef2f7"/>')
        body.append(svg_text(x, top + plot_h + 22, label, size=11, fill="#4b5563"))

    body.append(f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#111827"/>')
    body.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#111827"/>')
    body.append(svg_text(left + plot_w / 2, height - 28, xlabel, size=13, fill="#111827"))
    body.append(
        f'<text x="22" y="{top + plot_h / 2:.2f}" font-family="Arial, sans-serif" '
        f'font-size="13" text-anchor="middle" fill="#111827" '
        f'transform="rotate(-90 22 {top + plot_h / 2:.2f})">{html.escape(ylabel)}</text>'
    )

    for series in series_list:
        points = " ".join(f"{px(x):.2f},{py(y):.2f}" for x, y in series.points)
        body.append(f'<polyline points="{points}" fill="none" stroke="{series.color}" stroke-width="2.4"/>')
        for x, y in series.points:
            body.append(f'<circle cx="{px(x):.2f}" cy="{py(y):.2f}" r="4" fill="{series.color}" stroke="#ffffff" stroke-width="1"/>')

    legend_x = width - right + 32
    legend_y = top + 10
    for index, series in enumerate(series_list):
        y = legend_y + index * 22
        body.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x + 22}" y2="{y}" stroke="{series.color}" stroke-width="3"/>')
        body.append(f'<circle cx="{legend_x + 11}" cy="{y}" r="3.5" fill="{series.color}"/>')
        body.append(svg_text(legend_x + 30, y + 4, series.label, size=12, anchor="start", fill="#374151"))

    write_svg(path, width, height, body)


def selected_series(
    stats: list[GroupStats],
    value_name: str,
    transform_x,
    selected: list[tuple[str, str]],
) -> list[Series]:
    lookup = stat_lookup(stats)
    series_list: list[Series] = []
    sizes = sorted({item.L for item in stats})
    for model, k in selected:
        points = []
        for L in sizes:
            item = lookup.get((model, k, L))
            if item is None:
                continue
            points.append((transform_x(L), getattr(item, value_name)))
        if points:
            series_list.append(
                Series(
                    label=model_label(model, k),
                    points=points,
                    color=key_color(model, k),
                )
            )
    return series_list


def k_series(stats: list[GroupStats], model: str, value_name: str, transform_x) -> list[Series]:
    series_list: list[Series] = []
    sizes = sorted({item.L for item in stats})
    k_values = sorted(
        {item.k for item in stats if item.model == model and item.k != "None"},
        key=lambda item: float(item),
        reverse=True,
    )
    lookup = stat_lookup(stats)
    for k in k_values:
        points = [
            (transform_x(L), getattr(lookup[(model, k, L)], value_name))
            for L in sizes
            if (model, k, L) in lookup
        ]
        series_list.append(Series(label=f"k={fmt_k(k)}", points=points, color=key_color(model, k)))
    return series_list


def disorder_sweep(
    path: Path,
    stats: list[GroupStats],
    L: int,
    value_name: str = "var_winding",
) -> None:
    series_list = []
    for model in ("gamma", "gamma4"):
        points = []
        for item in stats:
            if item.model == model and item.L == L and item.k != "None":
                points.append((math.log(float(item.k)), getattr(item, value_name)))
        points.sort()
        if points:
            series_list.append(Series(model, points, COLORS[model]))
    x_ticks = [
        (math.log(k), f"{k:g}")
        for k in sorted({float(item.k) for item in stats if item.k != "None"})
    ]
    line_plot(
        path=path,
        title=f"Disorder Sweep At L={L}",
        xlabel="Gamma shape k (log axis; smaller k = stronger disorder)",
        ylabel="winding variance",
        series_list=series_list,
        x_ticks=x_ticks,
    )


def scatter_boundary_hits(
    path: Path,
    rows: list[dict[str, str]],
    selected: list[tuple[str, str]],
    L: int,
    max_points: int = 650,
) -> None:
    width, height = 960, 960
    margin = 40
    panel_gap = 24
    panel_w = (width - 2 * margin - panel_gap * (len(selected) - 1)) / len(selected)
    panel_h = height - 100
    panel_size = min(panel_w, panel_h)
    top = 58
    body: list[str] = [svg_text(width / 2, 30, f"Boundary Hit Locations, L={L}", size=19, weight="700")]

    rng = random.Random(123)
    for idx, (model, k) in enumerate(selected):
        x0 = margin + idx * (panel_w + panel_gap)
        y0 = top
        group = [
            row
            for row in rows
            if row["model"] == model and (row["k"] or "None") == k and int(row["L"]) == L
        ]
        if len(group) > max_points:
            group = rng.sample(group, max_points)

        body.append(f'<rect x="{x0:.2f}" y="{y0}" width="{panel_size:.2f}" height="{panel_size:.2f}" fill="#f9fafb" stroke="#d1d5db"/>')
        body.append(svg_text(x0 + panel_size / 2, y0 - 14, model_label(model, k), size=13, weight="700"))

        def px(value: float) -> float:
            return x0 + (value / L + 1) / 2 * panel_size

        def py(value: float) -> float:
            return y0 + (1 - (value / L + 1) / 2) * panel_size

        body.append(f'<line x1="{px(0):.2f}" y1="{y0}" x2="{px(0):.2f}" y2="{y0 + panel_size}" stroke="#e5e7eb"/>')
        body.append(f'<line x1="{x0:.2f}" y1="{py(0):.2f}" x2="{x0 + panel_size:.2f}" y2="{py(0):.2f}" stroke="#e5e7eb"/>')
        for row in group:
            x = int(row["hit_x"])
            y = int(row["hit_y"])
            body.append(f'<circle cx="{px(x):.2f}" cy="{py(y):.2f}" r="2.2" fill="{key_color(model, k)}" opacity="0.55"/>')

        mean_x = mean(int(row["hit_x"]) for row in group) if group else 0
        mean_y = mean(int(row["hit_y"]) for row in group) if group else 0
        body.append(f'<circle cx="{px(mean_x):.2f}" cy="{py(mean_y):.2f}" r="5.2" fill="#111827" stroke="#ffffff" stroke-width="1.2"/>')
        body.append(svg_text(x0 + panel_size / 2, y0 + panel_size + 24, f"mean=({mean_x / L:.2f}, {mean_y / L:.2f})", size=12, fill="#374151"))

    write_svg(path, width, height, body)


def winding_histograms(
    path: Path,
    rows: list[dict[str, str]],
    selected: list[tuple[str, str]],
    L: int,
) -> None:
    width, height = 960, 960
    margin = 42
    panel_gap = 24
    panel_w = (width - 2 * margin - panel_gap * (len(selected) - 1)) / len(selected)
    panel_h = height - 116
    top = 62
    body: list[str] = [svg_text(width / 2, 31, f"Winding Histograms, L={L}", size=19, weight="700")]

    groups: list[tuple[tuple[str, str], Counter[int]]] = []
    all_windings = []
    for model, k in selected:
        values = [
            int(row["winding"])
            for row in rows
            if row["model"] == model and (row["k"] or "None") == k and int(row["L"]) == L
        ]
        counter = Counter(values)
        groups.append(((model, k), counter))
        all_windings.extend(values)

    if not all_windings:
        return

    min_w, max_w = min(all_windings), max(all_windings)
    max_count = max((max(counter.values()) if counter else 0) for _, counter in groups)

    for idx, ((model, k), counter) in enumerate(groups):
        x0 = margin + idx * (panel_w + panel_gap)
        y0 = top
        body.append(f'<rect x="{x0:.2f}" y="{y0}" width="{panel_w:.2f}" height="{panel_h}" fill="#f9fafb" stroke="#d1d5db"/>')
        body.append(svg_text(x0 + panel_w / 2, y0 - 14, model_label(model, k), size=13, weight="700"))

        def px(value: float) -> float:
            return x0 + (value - min_w) / max(max_w - min_w, 1) * panel_w

        def py(value: float) -> float:
            return y0 + panel_h - value / max(max_count, 1) * panel_h

        bar_w = max(panel_w / max(max_w - min_w + 1, 1) * 0.75, 2)
        for winding_value in range(min_w, max_w + 1):
            count = counter.get(winding_value, 0)
            x = px(winding_value)
            y = py(count)
            body.append(
                f'<rect x="{x - bar_w / 2:.2f}" y="{y:.2f}" width="{bar_w:.2f}" '
                f'height="{y0 + panel_h - y:.2f}" fill="{key_color(model, k)}" opacity="0.75"/>'
            )
        body.append(svg_text(x0 + panel_w / 2, y0 + panel_h + 24, "winding", size=12, fill="#374151"))

    write_svg(path, width, height, body)


def panel_path_plot(path: Path, L: int, selected: list[tuple[str, float | None]], seed: int) -> None:
    width, height = 1200, 1200
    margin = 36
    columns = 3
    rows = math.ceil(len(selected) / columns)
    panel_gap_x = 30
    panel_gap_y = 58
    panel_w = (width - 2 * margin - panel_gap_x * (columns - 1)) / columns
    panel_h = (height - 130 - panel_gap_y * (rows - 1)) / rows
    top = 72
    body: list[str] = [svg_text(width / 2, 31, f"Sample Loop-Erased Paths, L={L}", size=19, weight="700")]

    for idx, (model, k) in enumerate(selected):
        row = idx // columns
        col = idx % columns
        x0 = margin + col * (panel_w + panel_gap_x)
        y0 = top + row * (panel_h + panel_gap_y)
        rng = random.Random(seed + idx * 1009)
        env = RandomEnvironment(rng, model=model, k=k)
        raw_path = random_walk_to_boundary(L, env)
        erased = loop_erase(raw_path)
        w = winding(erased)
        color = key_color(model, "None" if k is None else str(k))

        panel_size = min(panel_w, panel_h)
        body.append(f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{panel_size:.2f}" height="{panel_size:.2f}" fill="#f9fafb" stroke="#d1d5db"/>')
        body.append(svg_text(x0 + panel_size / 2, y0 - 10, model_label(model, None if k is None else str(k)), size=13, weight="700"))

        pad = 14
        def px(value: float) -> float:
            return x0 + pad + (value + L) / (2 * L) * (panel_size - 2 * pad)

        def py(value: float) -> float:
            return y0 + pad + (L - value) / (2 * L) * (panel_size - 2 * pad)

        for frac in (-0.5, 0.0, 0.5):
            gx = px(frac * L)
            gy = py(frac * L)
            body.append(f'<line x1="{gx:.2f}" y1="{y0 + pad:.2f}" x2="{gx:.2f}" y2="{y0 + panel_size - pad:.2f}" stroke="#e5e7eb"/>')
            body.append(f'<line x1="{x0 + pad:.2f}" y1="{gy:.2f}" x2="{x0 + panel_size - pad:.2f}" y2="{gy:.2f}" stroke="#e5e7eb"/>')

        points = " ".join(f"{px(x):.2f},{py(y):.2f}" for x, y in erased)
        body.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.1" stroke-linejoin="round" stroke-linecap="round"/>')
        sx, sy = erased[0]
        hx, hy = erased[-1]
        body.append(f'<circle cx="{px(sx):.2f}" cy="{py(sy):.2f}" r="4" fill="#111827"/>')
        body.append(f'<circle cx="{px(hx):.2f}" cy="{py(hy):.2f}" r="4" fill="{color}"/>')
        body.append(svg_text(x0 + panel_size / 2, y0 + panel_size + 18, f"W={w}, raw={len(raw_path)-1}, LERW={len(erased)-1}", size=11, fill="#374151"))

    write_svg(path, width, height, body)


def environment_vector_fields(path: Path, L: int, selected: list[tuple[str, float | None]], seed: int) -> None:
    width, height = 960, 960
    margin = 42
    panel_gap = 30
    panel_w = (width - 2 * margin - panel_gap * (len(selected) - 1)) / len(selected)
    panel_h = height - 100
    panel_size = min(panel_w, panel_h)
    top = 62
    body: list[str] = [svg_text(width / 2, 31, f"Local Mean-Step Fields, L={L}", size=19, weight="700")]

    for idx, (model, k) in enumerate(selected):
        x0 = margin + idx * (panel_w + panel_gap)
        y0 = top
        rng = random.Random(seed + idx * 2003)
        env = RandomEnvironment(rng, model=model, k=k)
        color = key_color(model, "None" if k is None else str(k))
        body.append(f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{panel_size:.2f}" height="{panel_size:.2f}" fill="#f9fafb" stroke="#d1d5db"/>')
        body.append(svg_text(x0 + panel_size / 2, y0 - 14, model_label(model, None if k is None else str(k)), size=13, weight="700"))

        pad = 16
        def px(value: float) -> float:
            return x0 + pad + (value + L) / (2 * L) * (panel_size - 2 * pad)

        def py(value: float) -> float:
            return y0 + pad + (L - value) / (2 * L) * (panel_size - 2 * pad)

        for x in range(-L + 1, L):
            for y in range(-L + 1, L):
                p_n, p_e, p_s, p_w = env.transition_probabilities((x, y))
                dx = p_e - p_w
                dy = p_n - p_s
                scale = min(panel_w, panel_h) * 0.12
                x1, y1 = px(x), py(y)
                x2, y2 = x1 + dx * scale, y1 - dy * scale
                body.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" stroke-width="1.15" opacity="0.65"/>')
                body.append(f'<circle cx="{x2:.2f}" cy="{y2:.2f}" r="1.3" fill="{color}" opacity="0.75"/>')

    write_svg(path, width, height, body)


def write_group_stats_csv(path: Path, stats: list[GroupStats]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "model",
                "k",
                "L",
                "samples",
                "mean_winding",
                "var_winding",
                "mean_raw_steps",
                "mean_lerw_steps",
                "mean_hit_x",
                "mean_hit_y",
                "mean_hit_x_over_L",
                "mean_hit_y_over_L",
            ]
        )
        for item in stats:
            writer.writerow(
                [
                    item.model,
                    item.k,
                    item.L,
                    item.samples,
                    item.mean_winding,
                    item.var_winding,
                    item.mean_raw_steps,
                    item.mean_lerw_steps,
                    item.mean_hit_x,
                    item.mean_hit_y,
                    item.mean_hit_x / item.L,
                    item.mean_hit_y / item.L,
                ]
            )


def write_fit_csv(path: Path, stats: list[GroupStats]) -> list[dict[str, str | float]]:
    by_model: dict[tuple[str, str], list[GroupStats]] = defaultdict(list)
    for item in stats:
        by_model[(item.model, item.k)].append(item)

    rows: list[dict[str, str | float]] = []
    for (model, k), items in sorted(by_model.items()):
        items = sorted(items, key=lambda item: item.L)
        if len(items) < 3:
            continue
        ys = [item.var_winding for item in items]
        for fit_name, xs in (
            ("logL", [math.log(item.L) for item in items]),
            ("logL_squared", [math.log(item.L) ** 2 for item in items]),
        ):
            slope, intercept, r2 = linear_fit(xs, ys)
            fitted = [slope * x + intercept for x in xs]
            rss = sum((y - y_hat) ** 2 for y, y_hat in zip(ys, fitted))
            rows.append(
                {
                    "model": model,
                    "k": k,
                    "fit": fit_name,
                    "slope": slope,
                    "intercept": intercept,
                    "r_squared": r2,
                    "rss": rss,
                }
            )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["model", "k", "fit", "slope", "intercept", "r_squared", "rss"],
        )
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_analysis_report(
    path: Path,
    source_csv: Path,
    stats: list[GroupStats],
    fit_rows: list[dict[str, str | float]],
    figure_names: list[str],
) -> None:
    largest_L = max(item.L for item in stats)
    stats_by_key = stat_lookup(stats)
    lines: list[str] = []
    lines.append("# Numerical Analysis Dossier")
    lines.append("")
    lines.append(f"Source data: `{source_csv}`")
    lines.append("")
    lines.append("## Dataset")
    lines.append("")
    models = sorted({item.model for item in stats})
    sizes = sorted({item.L for item in stats})
    lines.append(f"- Models: {', '.join(models)}")
    lines.append(f"- Box sizes: {', '.join(str(size) for size in sizes)}")
    lines.append(f"- Total grouped cells: {len(stats)}")
    lines.append("")
    lines.append("## Main Diagnostics")
    lines.append("")
    for model, k in [("symmetric", "None"), ("gamma", "1.0"), ("gamma4", "1.0")]:
        item = stats_by_key.get((model, k, largest_L))
        if item is None:
            continue
        lines.append(
            f"- `{model}` k={k} at L={largest_L}: "
            f"Var(W)={item.var_winding:.3f}, mean W={item.mean_winding:.3f}, "
            f"mean hit/L=({item.mean_hit_x / item.L:.3f}, {item.mean_hit_y / item.L:.3f}), "
            f"mean raw steps={item.mean_raw_steps:.1f}."
        )
    lines.append("")
    lines.append("The two-random-weight `gamma` model should be treated carefully: the boundary-hit diagnostics measure whether the normalised probabilities create an effective drift. The `gamma4` model is a balanced comparison, not a replacement for the supervisor-suggested model.")
    lines.append("")
    lines.append("## Scaling Fits")
    lines.append("")
    lines.append("The fits below are linear regressions of winding variance against either `log L` or `(log L)^2`. They are descriptive diagnostics, not proof of asymptotic scaling.")
    lines.append("")
    lines.append("| model | k | fit | slope | intercept | R^2 | RSS |")
    lines.append("|---|---:|---|---:|---:|---:|---:|")
    for row in fit_rows:
        lines.append(
            f"| {row['model']} | {row['k']} | {row['fit']} | "
            f"{float(row['slope']):.4g} | {float(row['intercept']):.4g} | "
            f"{float(row['r_squared']):.4f} | {float(row['rss']):.4g} |"
        )
    lines.append("")
    lines.append("## Figures")
    lines.append("")
    for name in figure_names:
        lines.append(f"![{name}](figures/{name})")
        lines.append("")
    lines.append("## Preliminary Interpretation")
    lines.append("")
    lines.append("- The symmetric baseline is the calibration case and should show roughly logarithmic growth, subject to finite-size noise.")
    lines.append("- A strong mean boundary-hit displacement indicates effective drift, in which case winding variance may flatten rather than grow like `log L` or `(log L)^2`.")
    lines.append("- The balanced `gamma4` diagnostic helps distinguish disorder effects from drift effects.")
    lines.append("- The current run should be used to decide which models and size ranges deserve larger sample counts before paper drafting.")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_all_figures(csv_path: Path, output_dir: Path) -> None:
    rows = read_rows([csv_path])
    stats = grouped_stats(rows)
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    write_group_stats_csv(tables_dir / "group_stats.csv", stats)
    fit_rows = write_fit_csv(tables_dir / "scaling_fits.csv", stats)

    sizes = sorted({item.L for item in stats})
    x_ticks_log = [(math.log(L), str(L)) for L in sizes]
    x_ticks_log2 = [(math.log(L) ** 2, str(L)) for L in sizes]
    selected = [
        ("symmetric", "None"),
        ("gamma", "1.0"),
        ("gamma4", "20.0"),
        ("gamma4", "1.0"),
        ("gamma4", "0.5"),
    ]

    figure_names: list[str] = []

    def add_line_plot(filename: str, **kwargs) -> None:
        line_plot(path=figures_dir / filename, **kwargs)
        figure_names.append(filename)

    add_line_plot(
        "variance_vs_logL_selected.svg",
        title="Winding Variance vs log L",
        xlabel="box size L (shown at log L positions)",
        ylabel="sample variance of winding",
        series_list=selected_series(stats, "var_winding", math.log, selected),
        x_ticks=x_ticks_log,
    )
    add_line_plot(
        "variance_vs_logL2_selected.svg",
        title="Winding Variance vs (log L)^2",
        xlabel="box size L (shown at (log L)^2 positions)",
        ylabel="sample variance of winding",
        series_list=selected_series(stats, "var_winding", lambda L: math.log(L) ** 2, selected),
        x_ticks=x_ticks_log2,
    )
    add_line_plot(
        "gamma_variance_by_k_logL.svg",
        title="Two-Weight Gamma Model: Variance By Disorder",
        xlabel="box size L (shown at log L positions)",
        ylabel="sample variance of winding",
        series_list=k_series(stats, "gamma", "var_winding", math.log),
        x_ticks=x_ticks_log,
    )
    add_line_plot(
        "gamma4_variance_by_k_logL.svg",
        title="Balanced Gamma4 Diagnostic: Variance By Disorder",
        xlabel="box size L (shown at log L positions)",
        ylabel="sample variance of winding",
        series_list=k_series(stats, "gamma4", "var_winding", math.log),
        x_ticks=x_ticks_log,
    )
    add_line_plot(
        "mean_hit_x_over_L.svg",
        title="Drift Diagnostic: Mean Boundary Hit x/L",
        xlabel="box size L (shown at log L positions)",
        ylabel="mean hit x / L",
        series_list=[
            Series(item.label, [(x, y / math.exp(x)) for x, y in item.points], item.color)
            for item in selected_series(stats, "mean_hit_x", math.log, selected)
        ],
        x_ticks=x_ticks_log,
    )
    add_line_plot(
        "mean_hit_y_over_L.svg",
        title="Drift Diagnostic: Mean Boundary Hit y/L",
        xlabel="box size L (shown at log L positions)",
        ylabel="mean hit y / L",
        series_list=[
            Series(item.label, [(x, y / math.exp(x)) for x, y in item.points], item.color)
            for item in selected_series(stats, "mean_hit_y", math.log, selected)
        ],
        x_ticks=x_ticks_log,
    )
    add_line_plot(
        "raw_walk_length.svg",
        title="Mean Raw Walk Length",
        xlabel="box size L (shown at log L positions)",
        ylabel="mean raw steps",
        series_list=selected_series(stats, "mean_raw_steps", math.log, selected),
        x_ticks=x_ticks_log,
    )
    add_line_plot(
        "lerw_path_length.svg",
        title="Mean Loop-Erased Path Length",
        xlabel="box size L (shown at log L positions)",
        ylabel="mean LERW steps",
        series_list=selected_series(stats, "mean_lerw_steps", math.log, selected),
        x_ticks=x_ticks_log,
    )

    largest_L = max(sizes)
    mid_L = 128 if 128 in sizes else largest_L
    disorder_sweep(figures_dir / f"disorder_sweep_L{mid_L}.svg", stats, mid_L)
    figure_names.append(f"disorder_sweep_L{mid_L}.svg")
    if largest_L != mid_L:
        disorder_sweep(figures_dir / f"disorder_sweep_L{largest_L}.svg", stats, largest_L)
        figure_names.append(f"disorder_sweep_L{largest_L}.svg")

    boundary_selected = [("symmetric", "None"), ("gamma", "1.0"), ("gamma4", "1.0")]
    scatter_boundary_hits(figures_dir / f"boundary_hits_L{mid_L}.svg", rows, boundary_selected, mid_L)
    figure_names.append(f"boundary_hits_L{mid_L}.svg")
    winding_histograms(figures_dir / f"winding_histograms_L{mid_L}.svg", rows, boundary_selected, mid_L)
    figure_names.append(f"winding_histograms_L{mid_L}.svg")
    panel_path_plot(
        figures_dir / "sample_lerw_paths.svg",
        L=64,
        selected=[
            ("symmetric", None),
            ("gamma", 1.0),
            ("gamma4", 1.0),
            ("gamma", 20.0),
            ("gamma4", 20.0),
            ("gamma4", 0.5),
        ],
        seed=20260612,
    )
    figure_names.append("sample_lerw_paths.svg")
    environment_vector_fields(
        figures_dir / "environment_vector_fields.svg",
        L=7,
        selected=[("symmetric", None), ("gamma", 1.0), ("gamma4", 1.0)],
        seed=20260612,
    )
    figure_names.append("environment_vector_fields.svg")

    write_analysis_report(
        output_dir / "analysis_dossier.md",
        source_csv=csv_path,
        stats=stats,
        fit_rows=fit_rows,
        figure_names=figure_names,
    )
    print(f"wrote {len(figure_names)} figures to {figures_dir}")
    print(f"wrote tables to {tables_dir}")
    print(f"wrote report to {output_dir / 'analysis_dossier.md'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate figures for LERW results.")
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("analysis_output"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    make_all_figures(args.csv_file, args.output_dir)


if __name__ == "__main__":
    main()
