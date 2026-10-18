import math
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pendulab as pl  # noqa: E402


def _period(length, g, amp):
    return pl.small_angle_period(length, g) * pl.amplitude_factor(amp)


def test_parse_valid_run_and_peaks():
    lab = pl.parse_lab(
        "# comment\n"
        "\n"
        "run len=0.5 n=20 amp=5 t=28.38,28.36 label=A\n"
        "peaks len=1.0 a=10,9,8 label=P\n"
    )
    assert lab.errors == []
    run = lab.runs[0]
    assert (run.line, run.length, run.amp, run.n, run.label) == (3, 0.5, 5.0, 20, "A")
    assert run.times == [28.38, 28.36]
    assert lab.peaks[0].amps == [10.0, 9.0, 8.0]


def test_parse_errors_carry_line_numbers():
    lab = pl.parse_lab(
        "run len=0.5 n=20 t=28.3\n"          # ok
        "run len=-1 t=10\n"                  # negative length
        "run len=1 t=abc\n"                  # bad number
        "run len=1 t=10 colour=red\n"        # unknown key
        "run len=1 n=2.5 t=10\n"             # n must be whole
        "run len=1 amp=95 t=10\n"            # amplitude out of range
        "run len=1\n"                        # missing t
        "bogus len=1\n"                      # unknown record
    )
    lines = sorted({line for line, _ in lab.errors})
    assert lines == [2, 3, 4, 5, 6, 7, 8]
    assert len(lab.runs) == 1


def test_peaks_must_strictly_decrease():
    lab = pl.parse_lab("peaks len=1 a=10,10.5,9\n")
    assert any("strictly decrease" in msg for _, msg in lab.errors)
    assert lab.peaks == []


def test_empty_file_reports_no_records():
    lab = pl.parse_lab("# nothing here\n")
    assert lab.errors == [(0, "file has no run or peaks records")]


def test_amplitude_factor_limits():
    assert pl.amplitude_factor(0.0) == pytest.approx(1.0, abs=1e-12)
    theta = math.radians(3.0)
    series = 1.0 + theta ** 2 / 16.0 + 11.0 * theta ** 4 / 3072.0
    assert pl.amplitude_factor(3.0) == pytest.approx(series, rel=1e-8)
    assert pl.ellipk(0.0) == pytest.approx(math.pi / 2, rel=1e-12)


def test_analyse_run_recovers_g_from_large_angle_swing():
    g_true = 9.81
    t = _period(0.75, g_true, 12.0)
    run = pl.Run(line=1, length=0.75, amp=12.0, n=1, times=[t, t, t], label="x")
    result = pl.analyse_run(run, dl=0.0)
    assert result.g == pytest.approx(g_true, rel=1e-9)
    assert result.sd_t == pytest.approx(0.0, abs=1e-12)


def test_chauvenet_rejects_gross_outlier_only():
    values = [10.0, 10.01, 9.99, 10.02, 9.98, 10.01, 9.99, 10.0, 12.0]
    keep = pl.chauvenet_keep(values)
    assert keep[-1] is False
    assert all(keep[:-1])


def test_gravity_formula_reference_values():
    assert pl.gravity_at(0.0) == pytest.approx(9.7803, abs=1e-4)
    assert pl.gravity_at(90.0) == pytest.approx(9.8322, abs=1e-4)
    assert pl.gravity_at(45.0, alt_m=1000.0) < pl.gravity_at(45.0)


def test_fit_through_origin_recovers_slope():
    g = 9.81
    lengths = [0.25, 0.5, 1.0, 1.5]
    k = 4 * math.pi ** 2 / g
    s, se, r2 = pl.fit_through_origin(lengths, [k * x for x in lengths])
    assert s == pytest.approx(k, rel=1e-12)
    assert r2 == pytest.approx(1.0)


def test_fit_line_recovers_intercept():
    xs = [0.25, 0.5, 1.0, 1.5]
    a, b, r2 = pl.fit_line(xs, [0.05 + 2.0 * x for x in xs])
    assert a == pytest.approx(0.05, abs=1e-12)
    assert b == pytest.approx(2.0, abs=1e-12)


def test_damping_from_geometric_peaks():
    amps = [10.0 * 0.9 ** i for i in range(11)]
    d = pl.damping(amps)
    assert d["cycles"] == 10
    assert d["delta"] == pytest.approx(-math.log(0.9), rel=1e-12)
    assert d["q"] == pytest.approx(math.pi / d["delta"])
    assert d["after_10"] == pytest.approx(0.9 ** 10, rel=1e-9)


def test_predict_and_length_round_trip():
    period = _period(0.9, 9.80, 8.0)
    # length() should invert predict() exactly for the same amplitude
    t0 = period / pl.amplitude_factor(8.0)
    assert 9.80 * (t0 / (2 * math.pi)) ** 2 == pytest.approx(0.9, rel=1e-12)


def _write(tmp_path, text):
    path = tmp_path / "lab.pend"
    path.write_text(text, encoding="utf-8")
    return str(path)


def test_cli_ci_gate_exit_codes(tmp_path, capsys):
    t = _period(1.0, 9.81, 5.0)
    path = _write(tmp_path, f"run len=1 amp=5 t={t!r},{t!r}\nrun len=0.5 amp=5 "
                            f"t={_period(0.5, 9.81, 5.0)!r},{_period(0.5, 9.81, 5.0)!r}\n")
    assert pl.main(["analyze", path, "--expect-g", "9.81", "--tol-pct", "0.1"]) == 0
    assert pl.main(["analyze", path, "--expect-g", "9.0", "--tol-pct", "1"]) == 2
    assert "FAIL" in capsys.readouterr().err


def test_cli_lint_fails_on_bad_file(tmp_path, capsys):
    path = _write(tmp_path, "run len=1\n")
    assert pl.main(["lint", path]) == 1
    out = capsys.readouterr().out
    assert f"{path}:1: error:" in out
