"""Mutation checks for scientific evidence drift, without recompiling LaTeX."""
import json
from pathlib import Path
import shutil

import pytest

import manuscript as m


@pytest.fixture
def copy(tmp_path):
    for path in [*m.input_paths(), "paper/build-report.json"]:
        target = tmp_path / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(m.ROOT / path, target)
    return tmp_path


def mutate(root, path, change):
    target = root / path
    data = json.loads(target.read_text())
    change(data)
    target.write_text(json.dumps(data))


def test_versioned_evidence_is_current():
    m.check_tables()


def test_reject_changed_numerical_source(copy):
    with (copy / "phasechange/material.py").open("a") as f:
        f.write("\n# changed constitutive implementation\n")
    with pytest.raises(ValueError, match="Stale numerical source"):
        m.check_tables(copy)


def test_reject_missing_source_hash(copy):
    mutate(copy, m.REPORTS[-1], lambda d: d["source_sha256"].pop("phasechange/continuum.py"))
    with pytest.raises(ValueError, match="Incomplete refinement source"):
        m.check_tables(copy)


def test_reject_new_untracked_numerical_source(copy):
    (copy / "phasechange/new_solver.py").write_text("# source inventory must be refreshed\n")
    with pytest.raises(ValueError, match="Incomplete refinement source"):
        m.check_tables(copy)


def test_reject_failed_gate(copy):
    mutate(copy, m.REPORTS[-1], lambda d: d["acceptance"].update(refined_stefan_front=False))
    with pytest.raises(ValueError, match="acceptance"):
        m.check_tables(copy)


@pytest.mark.parametrize("error", [None, float("nan"), -0.1])
def test_reject_invalid_front(copy, error):
    mutate(copy, m.REPORTS[-1], lambda d: d["cases"]["refined"]["frozen"]["samples"][0].update(relative_front_error=error))
    with pytest.raises(ValueError, match="front error"):
        m.check_tables(copy)


def test_reject_missing_observation(copy):
    mutate(copy, m.REPORTS[-1], lambda d: d["cases"]["refined"]["enthalpy"]["samples"].pop())
    with pytest.raises(ValueError, match="observation window"):
        m.check_tables(copy)


def test_reject_changed_measurement_with_old_table(copy):
    mutate(copy, m.REPORTS[1], lambda d: d["heat_cool_slab"][0].update(max_sampled_method_difference_K=1.0))
    with pytest.raises(ValueError, match="Stale manuscript table"):
        m.check_tables(copy)


def test_reject_edited_table(copy):
    (copy / m.GENERATED[0]).write_text("fabricated table")
    with pytest.raises(ValueError, match="Stale manuscript table"):
        m.check_tables(copy)


def test_reject_missing_build_input(copy):
    mutate(copy, "paper/build-report.json", lambda d: d["sha256"].pop("paper/paper.pdf"))
    with pytest.raises(ValueError, match="Incomplete manuscript input"):
        m.check_report(copy)


def test_reject_changed_pdf(copy):
    with (copy / "paper/paper.pdf").open("ab") as f:
        f.write(b"\n% changed artifact\n")
    with pytest.raises(ValueError, match="Stale manuscript build"):
        m.check_report(copy)


def test_reject_changed_page_count(copy):
    # Populate a fresh manifest so this test checks page-count consistency alone.
    report = {"pages": 6, "sha256": {p: m.digest(copy / p) for p in m.input_paths(copy)}}
    (copy / "paper/build-report.json").write_text(json.dumps(report))
    with pytest.raises(ValueError, match="page-count"):
        m.check_report(copy)
