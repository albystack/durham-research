import json

import pandas as pd

from square_glauber.audit import write_height_convention_audit


def test_height_convention_audit_writes_full_matrices_and_figures(tmp_path) -> None:
    report = write_height_convention_audit(tmp_path, sizes=(6, 8))
    assert report["colour_swap_effect"] == "global sign reversal"
    assert report["reference_shift_effect"] == "global additive constant"
    matrix = pd.read_csv(tmp_path / "height_matrices_long.csv")
    assert len(matrix) == 2 * ((6 - 1) ** 2 + (8 - 1) ** 2)
    for L in (6, 8):
        assert (tmp_path / f"L{L}_all_horizontal_height_matrix.csv").exists()
        assert (tmp_path / f"L{L}_all_vertical_height_matrix.csv").exists()
        assert (tmp_path / f"L{L}_extremal_height_audit.png").exists()
    metadata = json.loads((tmp_path / "height_convention_metadata.json").read_text())
    assert metadata["white_vertex_rule"] == "(row + column) % 2 == 0"
