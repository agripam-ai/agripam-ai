"""Excel and TSV input/output for BioBrick parts and the experiment tracker."""
from __future__ import annotations

import io
from pathlib import Path
from typing import Any

import pandas as pd

from . import biobricks

PART_COLUMNS = ["part_id", "name", "role", "sequence", "hosts", "source", "notes"]


def _clean(value: Any) -> str:
    if value is None or (isinstance(value, float) and value != value):
        return ""
    return str(value).strip()


def parts_from_frame(frame: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    """Validate a parts table. Returns (parts, issues); rows starting with EXAMPLE are ignored."""
    issues, parts, seen = [], [], set()
    frame = frame.copy()
    frame.columns = [_clean(c).lower().replace(" ", "_") for c in frame.columns]
    for col in ("part_id", "sequence"):
        if col not in frame.columns:
            issues.append({"level": "error", "where": "Parts", "message": f"The column '{col}' is missing."})
            return [], issues
    for i, row in frame.iterrows():
        pid = _clean(row.get("part_id"))
        if not pid or pid.upper().startswith("EXAMPLE"):
            continue
        if pid in seen:
            issues.append({"level": "error", "where": "Parts", "message": f"Part '{pid}' appears more than once."})
            continue
        seen.add(pid)
        role = _clean(row.get("role")).lower() or "other"
        if role not in biobricks.ROLES:
            issues.append({"level": "warning", "where": "Parts",
                           "message": f"{pid}: role '{role}' is not one of {', '.join(biobricks.ROLES)}; using 'other'."})
            role = "other"
        try:
            seq = biobricks.clean_sequence(_clean(row.get("sequence")))
        except ValueError as error:
            issues.append({"level": "error", "where": "Parts", "message": f"{pid}: {error}"})
            continue
        if not seq:
            issues.append({"level": "error", "where": "Parts", "message": f"{pid}: the sequence is empty."})
            continue
        parts.append({"part_id": pid, "name": _clean(row.get("name")) or pid, "role": role, "sequence": seq,
                      "hosts": _clean(row.get("hosts")), "source": _clean(row.get("source")) or "user-supplied",
                      "notes": _clean(row.get("notes"))})
    if not parts and not any(i["level"] == "error" for i in issues):
        issues.append({"level": "error", "where": "Parts", "message": "No parts found. Replace the EXAMPLE rows with your own."})
    return parts, issues


def read_parts_workbook(source) -> tuple[list[dict], list[dict]]:
    sheets = pd.read_excel(source, sheet_name=None, dtype=object, engine="openpyxl")
    name = next((n for n in sheets if n.strip().lower() == "parts"), None)
    if name is None:
        return [], [{"level": "error", "where": "Parts", "message": "The sheet 'Parts' is missing."}]
    return parts_from_frame(sheets[name])


def read_parts_tsv(path: Path) -> list[dict]:
    parts, _ = parts_from_frame(pd.read_csv(path, sep="\t", dtype=object))
    return parts


def tracker_to_xlsx(rows: list[dict], analysis_notes: list[str] | None = None) -> bytes:
    buffer = io.BytesIO()
    frame = pd.DataFrame(rows, columns=biobricks.TRACKER_COLUMNS)
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        frame.to_excel(writer, sheet_name="Tracker", index=False)
        pd.DataFrame({"how to use": [
            "One row per construct or control. Fill the stage_* columns with a date when a stage is done.",
            "Record assay, replicates (3 or more), measured_value with units, and normalized_effect (test / parental).",
            "Every group needs a parental and a negative control with measured values.",
            "Upload this file again in the app to see progress and predicted-vs-measured agreement.",
        ] + (analysis_notes or [])}).to_excel(writer, sheet_name="Instructions", index=False)
    return buffer.getvalue()


def read_tracker(source) -> list[dict]:
    frame = pd.read_excel(source, sheet_name="Tracker", dtype=object, engine="openpyxl")
    frame = frame.where(pd.notna(frame), "")
    for col in biobricks.TRACKER_COLUMNS:
        if col not in frame.columns:
            frame[col] = ""
    return frame[biobricks.TRACKER_COLUMNS].to_dict(orient="records")
