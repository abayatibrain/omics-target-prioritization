"""Smoke tests for the Typer CLI."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from omics_target_prioritization.cli import app

runner = CliRunner()


def test_version() -> None:
    result = runner.invoke(app, ["version"])
    assert result.exit_code == 0
    assert "omics-target-prioritization" in result.stdout


def test_run_all_writes_outputs(tmp_path: Path) -> None:
    out = tmp_path / "out"
    result = runner.invoke(app, ["run-all", "--out", str(out), "--seed", "5"])
    assert result.exit_code == 0, result.stdout
    assert (out / "locus.json").exists()
    assert (out / "scores.json").exists()
    assert (out / "dossier.html").exists()


def test_step_by_step(tmp_path: Path) -> None:
    locus = tmp_path / "locus.json"
    evidence = tmp_path / "evidence.json"
    scores = tmp_path / "scores.json"
    dossier = tmp_path / "dossier.html"

    assert runner.invoke(app, ["simulate", "--out", str(locus)]).exit_code == 0
    assert (
        runner.invoke(app, ["integrate", "--locus", str(locus), "--out", str(evidence)]).exit_code
        == 0
    )
    assert (
        runner.invoke(app, ["score", "--evidence", str(evidence), "--out", str(scores)]).exit_code
        == 0
    )
    assert (
        runner.invoke(app, ["report", "--scores", str(scores), "--out", str(dossier)]).exit_code
        == 0
    )
    assert dossier.exists()
