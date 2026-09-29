import csv
from pathlib import Path


def test_confirmed_operons_are_read_from_confirmed_table(tmp_path: Path):
    confirmed = tmp_path / "cas_operons.tab"
    rejected = tmp_path / "cas_operons_putative.tab"
    confirmed.write_text("Prediction\tComplete_Interference\tBest_type\nI-C\t100%\tI-C\n")
    rejected.write_text("Prediction\tComplete_Interference\tBest_type\nFalse\t0%\tI-A\n")
    with confirmed.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["Best_type"] == "I-C"
    assert rows[0]["Complete_Interference"] == "100%"
