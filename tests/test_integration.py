from integration import experiment

def test_shared_model_contract():
    result = experiment()
    assert result["mass_operator_max_difference"] < 1e-12
    assert result["conductance_max_difference"] < 1e-12
    assert result["modal_reference_max_difference"] < 1e-10
    assert max(row["heat_balance_error"] for row in result["rows"]) < 1e-10
    for method in ("be", "cn", "rannacher"):
        rows = [r for r in result["rows"] if r["method"] == method]
        assert rows[-1]["mass_weighted_error"] < rows[0]["mass_weighted_error"]
