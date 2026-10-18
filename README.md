# Pendulab

**Offline simple-pendulum lab toolkit: timing runs in, gravity, fits, damping and CI gates out.**

Pendulab reads your stopwatch timings from one plain-text file and works out
the local value of *g*, with honest uncertainty. It applies the exact
large-angle correction, fits T² against length, measures damping from peak
amplitudes, and can gate a lab report in CI. It is written in pure Python
standard library: no network, no accounts, no API keys.

## Features

- **Plain-text lab files**: one `run` line per timing set, one `peaks` line per damping series, `#` comments.
- **Line-numbered errors**: every problem reports `file:line: error: message`, including unknown keys, bad numbers, out-of-range amplitudes and missing fields.
- **Lint command** with `--strict` to fail on warnings too, for pre-commit hooks.
- **Per-run g** from the period `T = total / n`, with sample standard deviation and standard error of the mean.
- **Exact large-angle correction** using the complete elliptic integral K (computed by the arithmetic-geometric mean), not just a series.
- **Outlier rejection** with Chauvenet's criterion; rejected timings are listed by line and index.
- **Uncertainty propagation** for g from length uncertainty (`--dl`) and timing scatter.
- **Overall g** from a least-squares fit of T₀² = k·L through the origin, with standard error and R².
- **Free-intercept fit** of T₀² = a + b·L, so a pivot offset shows up as a non-zero intercept.
- **Damping analysis**: logarithmic decrement, per-cycle spread, damping ratio ζ, Q, half-life in cycles and amplitude left after ten cycles.
- **Predictions**: period for a given length and amplitude, and the length needed for a target period.
- **Reference gravity** from latitude and altitude (1980 normal-gravity formula with free-air correction).
- **CSV export** of per-run results for spreadsheets or a lab notebook (`--csv`).
- **CI gate**: `analyze --expect-g 9.81 --tol-pct 0.5` exits with code 2 when the measured *g* is out of tolerance.

## Quick start

Requires Python 3.9 or later. No install step.

```bash
python3 pendulab.py analyze examples/bench.pend
```

`examples/bench.pend` holds four timing sets (0.25 m to 1.5 m, 20 swings each)
and one damping series. The `analyze` command prints a table of T, T₀ and *g*
per run, then the overall *g* with its uncertainty and the difference from
standard gravity.

Write the input yourself:

```text
# my_lab.pend
run len=0.750 n=20 amp=8 t=34.12,34.09,34.15 label=steel-bob
run len=1.200 n=20 amp=8 t=44.02,44.05,43.98 label=brass-bob
peaks len=1.200 a=8,7.7,7.4,7.1 label=damp
```

Then:

```bash
python3 pendulab.py analyze my_lab.pend --csv results.csv
python3 pendulab.py fit my_lab.pend
python3 pendulab.py damp my_lab.pend
```

## How to run

All commands take a file path (or numeric options) and print to stdout.

| Command | What it does |
| --- | --- |
| `analyze FILE [--csv PATH] [--dl M] [--expect-g G --tol-pct P]` | Per-run and overall *g*; optional CSV and CI gate |
| `fit FILE [--dl M]` | Origin fit and free-intercept fit of T₀² against length |
| `damp FILE` | Damping ratio, Q and half-life for each `peaks` record |
| `lint FILE [--strict]` | Check the file only; exit 1 on errors |
| `predict --len L [--amp A] [--n-osc N] (--g G \| --lat φ [--alt h])` | Period and stopwatch time for *N* swings |
| `length --period T [--amp A] (--g G \| --lat φ [--alt h])` | Length that gives a target period |
| `gravity --lat φ [--alt h]` | Normal gravity for a site |

Exit codes: `0` success, `1` input or usage error, `2` the `--expect-g` gate failed.

Without `--g` or `--lat`, predictions use standard gravity (9.80665 m/s²).

### File format

| Record | Keys | Meaning |
| --- | --- | --- |
| `run` | `len` (m, required), `t` (required, comma list of seconds), `n` (default 1), `amp` (degrees, default 5), `label` | Each `t` value is the total time for `n` swings |
| `peaks` | `len` (m, required), `a` (required, comma list of degrees, strictly decreasing), `label` | Successive peak amplitudes for damping |

Blank lines and lines starting with `#` are ignored. `label` is a single
token with no spaces.

## How to test

Tests use pytest and cover the parser, the physics, the fits, the damping
figures and the CLI exit codes.

```bash
python3 -m pip install pytest
python3 -m pytest tests -q
```

## Method notes

- Period per swing: `T = total / n`. The mean of the kept timings gives the measured period.
- Large-angle correction: `T = T₀ · 2K(sin(θ₀/2))/π`, so `T₀ = T / factor`, and `g = 4π²L / T₀²`.
- Uncertainty: `δg/g = √((δL/L)² + (2·sem/T)²)`. The amplitude uncertainty is not included.
- Origin fit: `k = Σ(L·T₀²) / Σ(L²)`, `g = 4π²/k`, with the standard error from the residuals (m−1 degrees of freedom).
- Damping: `δ` is the mean of `ln(aᵢ/aᵢ₊₁)`, `ζ = δ / √(4π² + δ²)`, `Q = π/δ`.
