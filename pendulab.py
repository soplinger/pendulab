#!/usr/bin/env python3
"""Pendulab: offline simple-pendulum lab toolkit (standard library only).

Plain-text timing runs in; gravity estimates, least-squares fits, damping
figures, predictions and CI gates out. No network, no third-party packages.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from dataclasses import dataclass, field

G_STANDARD = 9.80665     # standard gravity, m/s²
DEFAULT_DL = 0.001       # length uncertainty in metres (tape or ruler reading)
CHAUVENET_CUTOFF = 0.5   # reject a timing when N · P(|z| ≥ z_i) < 0.5
LARGE_AMP_DEG = 15.0     # above this, the large-angle correction matters
MAX_AMP_DEG = 90.0       # the simple-pendulum model is only offered below this
SPREAD_WARN_PCT = 1.0    # warn when timings spread more than this (CV, %)

RUN_KEYS = {"len", "n", "amp", "t", "label"}
PEAK_KEYS = {"len", "a", "label"}


# ---------------------------------------------------------------------------
# Input format
# ---------------------------------------------------------------------------
# One record per line. Blank lines and lines starting with '#' are ignored.
#   run   len=0.500 n=20 amp=5 t=28.38,28.36,28.41 label=bench-A
#   peaks len=0.500 a=10,9.7,9.4,9.1 label=damp-1
# 'len' is metres, 'amp' the release angle in degrees, 'n' the number of
# oscillations counted in each timing, and 't' the total seconds per timing.

@dataclass
class Run:
    line: int
    length: float
    amp: float
    n: int
    times: list
    label: str


@dataclass
class Peaks:
    line: int
    length: float
    amps: list
    label: str


@dataclass
class Lab:
    runs: list = field(default_factory=list)
    peaks: list = field(default_factory=list)
    errors: list = field(default_factory=list)     # (line, message)
    warnings: list = field(default_factory=list)   # (line, message)


def _pairs(tokens, line_no, allowed, errors):
    kv = {}
    for tok in tokens:
        key, sep, value = tok.partition("=")
        if not sep or not key or not value:
            errors.append((line_no, f"expected key=value, got '{tok}'"))
        elif key not in allowed:
            errors.append((line_no, f"unknown key '{key}'"))
        elif key in kv:
            errors.append((line_no, f"duplicate key '{key}'"))
        else:
            kv[key] = value
    return kv


def _number(text, name, line_no, errors):
    try:
        value = float(text)
    except ValueError:
        value = math.nan
    if not math.isfinite(value):
        errors.append((line_no, f"'{name}' is not a number: '{text}'"))
        return None
    return value


def _parse_run(line_no, tokens, lab):
    errors = lab.errors
    before = len(errors)
    kv = _pairs(tokens, line_no, RUN_KEYS, errors)
    for key in ("len", "t"):
        if key not in kv:
            errors.append((line_no, f"run is missing '{key}='"))
    if len(errors) > before:
        return

    length = _number(kv["len"], "len", line_no, errors)
    if length is not None and length <= 0:
        errors.append((line_no, f"'len' must be positive, got {kv['len']}"))

    amp = 5.0
    if "amp" in kv:
        amp = _number(kv["amp"], "amp", line_no, errors)
        if amp is not None and not 0 < amp < MAX_AMP_DEG:
            errors.append((line_no, f"'amp' must be between 0 and {MAX_AMP_DEG:g} degrees"))

    n_text = kv.get("n", "1")
    n = int(n_text) if n_text.isdigit() else 0
    if n < 1:
        errors.append((line_no, f"'n' must be a whole number of oscillations ≥ 1, got '{n_text}'"))

    times = []
    for part in kv["t"].split(","):
        value = _number(part, "t", line_no, errors)
        if value is not None and value <= 0:
            errors.append((line_no, f"timing must be positive, got {part}"))
        elif value is not None:
            times.append(value)

    if len(errors) > before:
        return
    if len(times) == 1:
        lab.warnings.append((line_no, "only one timing: no spread estimate for this run"))
    if amp > LARGE_AMP_DEG:
        lab.warnings.append((line_no, f"amplitude {amp:g}° is large: exact large-angle correction applied"))
    lab.runs.append(Run(line_no, length, amp, n, times, kv.get("label", "")))


def _parse_peaks(line_no, tokens, lab):
    errors = lab.errors
    before = len(errors)
    kv = _pairs(tokens, line_no, PEAK_KEYS, errors)
    for key in ("len", "a"):
        if key not in kv:
            errors.append((line_no, f"peaks record is missing '{key}='"))
    if len(errors) > before:
        return

    length = _number(kv["len"], "len", line_no, errors)
    if length is not None and length <= 0:
        errors.append((line_no, f"'len' must be positive, got {kv['len']}"))

    amps = []
    for part in kv["a"].split(","):
        value = _number(part, "a", line_no, errors)
        if value is not None and value <= 0:
            errors.append((line_no, f"peak amplitude must be positive, got {part}"))
        elif value is not None:
            amps.append(value)

    if len(amps) < 2 and len(errors) == before:
        errors.append((line_no, "peaks needs at least two amplitudes"))
    for i in range(1, len(amps)):
        if amps[i] >= amps[i - 1]:
            errors.append((line_no, "peak amplitudes must strictly decrease"))
            break

    if len(errors) > before:
        return
    lab.peaks.append(Peaks(line_no, length, amps, kv.get("label", "")))


def parse_lab(text):
    """Parse a pendulab file into runs, peak records, errors and warnings."""
    lab = Lab()
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        kind, *tokens = line.split()
        if kind == "run":
            _parse_run(line_no, tokens, lab)
        elif kind == "peaks":
            _parse_peaks(line_no, tokens, lab)
        else:
            lab.errors.append((line_no, f"unknown record '{kind}' (use 'run' or 'peaks')"))
    if not lab.runs and not lab.peaks and not lab.errors:
        lab.errors.append((0, "file has no run or peaks records"))
    return lab


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------

def ellipk(k):
    """Complete elliptic integral of the first kind K(k), by the AGM."""
    a, b = 1.0, math.sqrt(1.0 - k * k)
    while abs(a - b) > 1e-15 * a:
        a, b = (a + b) / 2.0, math.sqrt(a * b)
    return math.pi / (2.0 * a)


def amplitude_factor(amp_deg):
    """Exact ratio T / T0 for a swing of amp_deg (1.0 for small angles)."""
    k = math.sin(math.radians(amp_deg) / 2.0)
    return 2.0 * ellipk(k) / math.pi


def small_angle_period(length, g):
    return 2.0 * math.pi * math.sqrt(length / g)


def gravity_at(lat_deg, alt_m=0.0):
    """Normal gravity on the reference ellipsoid (1980 formula) plus free-air."""
    s1 = math.sin(math.radians(lat_deg)) ** 2
    s2 = math.sin(math.radians(2.0 * lat_deg)) ** 2
    g0 = 9.780327 * (1.0 + 0.0053024 * s1 - 0.0000058 * s2)
    return g0 - 3.086e-6 * alt_m


def chauvenet_keep(values):
    """Chauvenet's criterion: True for each value that is kept."""
    n = len(values)
    if n < 3:
        return [True] * n
    mean = statistics.fmean(values)
    sd = statistics.stdev(values)
    if sd == 0:
        return [True] * n
    keep = []
    for v in values:
        z = abs(v - mean) / sd
        keep.append(n * math.erfc(z / math.sqrt(2.0)) >= CHAUVENET_CUTOFF)
    return keep


@dataclass
class RunResult:
    line: int
    label: str
    length: float
    amp: float
    n: int
    timings: int
    used: int
    flagged: list
    mean_t: float
    sd_t: float | None
    t0: float
    g: float
    g_sd: float | None


def analyse_run(run, dl):
    per = [t / run.n for t in run.times]
    keep = chauvenet_keep(per)
    used = [p for p, k in zip(per, keep) if k] or per
    mean_t = statistics.fmean(used)
    sd_t = statistics.stdev(used) if len(used) >= 2 else None
    t0 = mean_t / amplitude_factor(run.amp)
    g = 4.0 * math.pi ** 2 * run.length / t0 ** 2
    g_sd = None
    if sd_t is not None:
        sem = sd_t / math.sqrt(len(used))
        rel = math.sqrt((dl / run.length) ** 2 + (2.0 * sem / mean_t) ** 2)
        g_sd = g * rel
    flagged = [i + 1 for i, k in enumerate(keep) if not k]
    return RunResult(run.line, run.label, run.length, run.amp, run.n, len(per),
                     len(used), flagged, mean_t, sd_t, t0, g, g_sd)


def fit_through_origin(xs, ys):
    """Least-squares slope of y = s·x. Returns (s, standard error, R²)."""
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    s = sxy / sxx
    ss_res = sum((y - s * x) ** 2 for x, y in zip(xs, ys))
    mean_y = statistics.fmean(ys)
    ss_tot = sum((y - mean_y) ** 2 for y in ys)
    se = math.sqrt(ss_res / (len(xs) - 1) / sxx) if len(xs) >= 2 else None
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return s, se, r2


def fit_line(xs, ys):
    """Ordinary least squares y = a + b·x. Returns (a, b, R²)."""
    m = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    denom = m * sxx - sx * sx
    b = (m * sxy - sx * sy) / denom
    a = (sy - b * sx) / m
    mean_y = sy / m
    ss_res = sum((y - a - b * x) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - mean_y) ** 2 for y in ys)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return a, b, r2


def damping(amps):
    """Logarithmic-decrement damping figures from successive peak amplitudes."""
    cycles = len(amps) - 1
    per_cycle = [math.log(amps[i] / amps[i + 1]) for i in range(cycles)]
    delta = sum(per_cycle) / cycles
    zeta = delta / math.sqrt(4.0 * math.pi ** 2 + delta ** 2)
    return {
        "cycles": cycles,
        "delta": delta,
        "min": min(per_cycle),
        "max": max(per_cycle),
        "zeta": zeta,
        "q": math.pi / delta,
        "half_life": math.log(2.0) / delta,
        "after_10": math.exp(-10.0 * delta),
    }


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _fmt(value, spec=".4f"):
    return "-" if value is None else format(value, spec)


def _render(headers, rows):
    cols = list(zip(headers, *rows)) if rows else [(h,) for h in headers]
    widths = [max(len(str(c)) for c in col) for col in cols]
    out = ["  ".join(str(h).ljust(w) for h, w in zip(headers, widths))]
    for row in rows:
        out.append("  ".join(str(c).ljust(w) for c, w in zip(row, widths)))
    return "\n".join(out)


def _diagnostics(path, items, kind):
    return [f"{path}:{line}: {kind}: {msg}" for line, msg in items]


def _overall_g(results):
    """Origin-fit g across distinct lengths, else the mean of the run values."""
    lengths = {round(r.length, 9) for r in results}
    if len(lengths) >= 2:
        s, se, r2 = fit_through_origin([r.length for r in results],
                                       [r.t0 ** 2 for r in results])
        g = 4.0 * math.pi ** 2 / s
        g_sd = g * se / s if se is not None else None
        return g, g_sd, "origin fit of T0² = k·L"
    gs = [r.g for r in results]
    g = statistics.fmean(gs)
    g_sd = statistics.stdev(gs) / math.sqrt(len(gs)) if len(gs) >= 2 else None
    return g, g_sd, "mean of runs (needs two lengths for a fit)"


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def _load(path):
    with open(path, encoding="utf-8") as fh:
        return parse_lab(fh.read())


def _report_parse(path, lab):
    for line in _diagnostics(path, sorted(lab.errors), "error"):
        print(line, file=sys.stderr)
    return 1 if lab.errors else 0


def cmd_lint(args):
    lab = _load(args.file)
    for line in _diagnostics(args.file, sorted(lab.errors), "error"):
        print(line)
    for line in _diagnostics(args.file, sorted(lab.warnings), "warning"):
        print(line)
    if lab.errors:
        print(f"{len(lab.errors)} error(s), {len(lab.warnings)} warning(s)")
        return 1
    print(f"ok: {len(lab.runs)} run(s), {len(lab.peaks)} peaks record(s), "
          f"{len(lab.warnings)} warning(s)")
    if args.strict and lab.warnings:
        return 1
    return 0


def cmd_analyze(args):
    lab = _load(args.file)
    if _report_parse(args.file, lab):
        return 1
    if not lab.runs:
        print(f"{args.file}: no run records to analyse", file=sys.stderr)
        return 1

    results = [analyse_run(r, args.dl) for r in lab.runs]
    g, g_sd, method = _overall_g(results)

    rows = []
    for r in results:
        rows.append([
            r.label or f"line {r.line}",
            f"{r.length:.4f}", f"{r.amp:g}", r.n,
            f"{r.used}/{r.timings}",
            _fmt(r.mean_t), _fmt(r.sd_t), _fmt(r.t0),
            _fmt(r.g, ".4f"), _fmt(r.g_sd, ".4f"),
        ])
    print(f"Pendulab analysis of {args.file}")
    print(_render(["run", "L (m)", "amp (°)", "n", "kept", "T (s)", "sd (s)",
                   "T0 (s)", "g (m/s²)", "±g"], rows))
    print()
    print(f"g = {g:.4f} ± {_fmt(g_sd)} m/s²   ({method})")
    print(f"Reference standard gravity {G_STANDARD} m/s²: "
          f"difference {g - G_STANDARD:+.4f} m/s² ({(g - G_STANDARD) / G_STANDARD * 100:+.3f} %)")

    flagged = [(r, i) for r in results for i in r.flagged]
    for r, i in flagged:
        print(f"note: line {r.line} timing {i} rejected by Chauvenet's criterion")
    for r in results:
        if r.sd_t is not None and r.mean_t and r.sd_t / r.mean_t * 100 > SPREAD_WARN_PCT:
            print(f"warning: line {r.line} spread is {r.sd_t / r.mean_t * 100:.2f} % "
                  f"of the mean period; recheck the timing")

    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["label", "length_m", "amp_deg", "n", "kept", "timings",
                             "mean_period_s", "sd_period_s", "T0_s", "g_m_s2", "g_sd_m_s2"])
            for r in results:
                writer.writerow([r.label, r.length, r.amp, r.n, r.used, r.timings,
                                 r.mean_t, r.sd_t if r.sd_t is not None else "",
                                 r.t0, r.g, r.g_sd if r.g_sd is not None else ""])
        print(f"wrote {args.csv}")

    if args.expect_g is not None:
        diff_pct = abs(g - args.expect_g) / args.expect_g * 100.0
        if diff_pct > args.tol_pct:
            print(f"FAIL: g {g:.4f} differs from expected {args.expect_g:g} "
                  f"by {diff_pct:.3f} % (limit {args.tol_pct:g} %)", file=sys.stderr)
            return 2
        print(f"PASS: g within {args.tol_pct:g} % of {args.expect_g:g}")
    return 0


def cmd_fit(args):
    lab = _load(args.file)
    if _report_parse(args.file, lab):
        return 1
    results = [analyse_run(r, args.dl) for r in lab.runs]
    xs = [r.length for r in results]
    ys = [r.t0 ** 2 for r in results]
    if len({round(x, 9) for x in xs}) < 2:
        print("fit needs runs at two or more different lengths", file=sys.stderr)
        return 1

    s, se, r2 = fit_through_origin(xs, ys)
    g = 4.0 * math.pi ** 2 / s
    g_sd = g * se / s if se is not None else None
    print(f"Fit of T0² = k·L through the origin ({len(xs)} runs)")
    print(f"  k = {s:.6f} s²/m   g = {g:.4f} ± {_fmt(g_sd)} m/s²   R² = {r2:.6f}")

    a, b, r2_line = fit_line(xs, ys)
    g_line = 4.0 * math.pi ** 2 / b
    print("Fit of T0² = a + b·L (free intercept)")
    print(f"  b = {b:.6f} s²/m   g = {g_line:.4f} m/s²   "
          f"intercept a = {a:+.6f} s²   R² = {r2_line:.6f}")
    print("  A large intercept points to a pivot offset or an unmeasured bob radius.")
    return 0


def cmd_damp(args):
    lab = _load(args.file)
    if _report_parse(args.file, lab):
        return 1
    if not lab.peaks:
        print(f"{args.file}: no peaks records to analyse", file=sys.stderr)
        return 1
    for rec in lab.peaks:
        d = damping(rec.amps)
        name = rec.label or f"line {rec.line}"
        print(f"{name}  L = {rec.length:g} m  cycles = {d['cycles']}")
        print(f"  log decrement δ = {d['delta']:.5f} "
              f"(per cycle {d['min']:.5f} to {d['max']:.5f})")
        print(f"  damping ratio ζ = {d['zeta']:.5f}   Q = {d['q']:.2f}")
        print(f"  amplitude half-life = {d['half_life']:.2f} cycles   "
              f"left after 10 cycles = {d['after_10'] * 100:.1f} %")
    return 0


def _resolve_g(args, parser):
    if args.g is not None and args.lat is not None:
        parser.error("use either --g or --lat, not both")
    if args.lat is not None:
        return gravity_at(args.lat, args.alt)
    return args.g if args.g is not None else G_STANDARD


def cmd_predict(args, parser):
    if args.len <= 0:
        parser.error("--len must be positive")
    if not 0 < args.amp < MAX_AMP_DEG:
        parser.error(f"--amp must be between 0 and {MAX_AMP_DEG:g} degrees")
    g = _resolve_g(args, parser)
    t = small_angle_period(args.len, g) * amplitude_factor(args.amp)
    print(f"L = {args.len:g} m   amp = {args.amp:g}°   g = {g:.5f} m/s²")
    print(f"period T = {t:.5f} s   (small-angle {small_angle_period(args.len, g):.5f} s)")
    print(f"{args.n_osc} oscillations take {t * args.n_osc:.3f} s")
    return 0


def cmd_length(args, parser):
    if args.period <= 0:
        parser.error("--period must be positive")
    if not 0 < args.amp < MAX_AMP_DEG:
        parser.error(f"--amp must be between 0 and {MAX_AMP_DEG:g} degrees")
    g = _resolve_g(args, parser)
    t0 = args.period / amplitude_factor(args.amp)
    length = g * (t0 / (2.0 * math.pi)) ** 2
    print(f"target T = {args.period:g} s at amp {args.amp:g}°, g = {g:.5f} m/s²")
    print(f"length L = {length:.5f} m  ({length * 1000:.2f} mm)")
    return 0


def cmd_gravity(args):
    g = gravity_at(args.lat, args.alt)
    print(f"lat {args.lat:g}°  alt {args.alt:g} m  ->  g = {g:.5f} m/s²")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        prog="pendulab",
        description="Offline simple-pendulum lab toolkit: timings in, g and damping out.")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("analyze", help="per-run and overall g from timing runs")
    p.add_argument("file")
    p.add_argument("--dl", type=float, default=DEFAULT_DL,
                   help="length uncertainty in metres (default 0.001)")
    p.add_argument("--csv", metavar="PATH", help="write per-run results as CSV")
    p.add_argument("--expect-g", type=float, help="CI gate: expected g in m/s²")
    p.add_argument("--tol-pct", type=float, default=1.0,
                   help="allowed deviation from --expect-g in percent (default 1)")
    p.set_defaults(func=cmd_analyze)

    p = sub.add_parser("fit", help="least-squares T0² against length")
    p.add_argument("file")
    p.add_argument("--dl", type=float, default=DEFAULT_DL)
    p.set_defaults(func=cmd_fit)

    p = sub.add_parser("damp", help="damping ratio and Q from peak amplitudes")
    p.add_argument("file")
    p.set_defaults(func=cmd_damp)

    p = sub.add_parser("lint", help="check a file and report line-numbered problems")
    p.add_argument("file")
    p.add_argument("--strict", action="store_true", help="treat warnings as failures")
    p.set_defaults(func=cmd_lint)

    p = sub.add_parser("predict", help="period for a given length and amplitude")
    p.add_argument("--len", type=float, required=True, help="length in metres")
    p.add_argument("--amp", type=float, default=5.0, help="amplitude in degrees")
    p.add_argument("--n-osc", type=int, default=20, help="oscillations for the stopwatch time")
    _add_gravity_options(p)
    p.set_defaults(func=cmd_predict, wants_parser=True)

    p = sub.add_parser("length", help="length that gives a target period")
    p.add_argument("--period", type=float, required=True, help="target period in seconds")
    p.add_argument("--amp", type=float, default=5.0, help="amplitude in degrees")
    _add_gravity_options(p)
    p.set_defaults(func=cmd_length, wants_parser=True)

    p = sub.add_parser("gravity", help="normal gravity for a latitude and altitude")
    p.add_argument("--lat", type=float, required=True, help="latitude in degrees")
    p.add_argument("--alt", type=float, default=0.0, help="altitude in metres")
    p.set_defaults(func=cmd_gravity)
    return parser


def _add_gravity_options(p):
    p.add_argument("--g", type=float, help="gravity in m/s² (default 9.80665)")
    p.add_argument("--lat", type=float, help="latitude in degrees, to derive g")
    p.add_argument("--alt", type=float, default=0.0, help="altitude in metres (with --lat)")


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "wants_parser", False):
        return args.func(args, parser)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
