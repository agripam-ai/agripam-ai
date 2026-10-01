"""State helpers for the session-scoped agricultural editing knowledgebase."""

from __future__ import annotations

import hashlib

import pandas as pd


RECORD_KEY_COLUMNS = ("organism", "strain", "editing_route", "source_id")


def normalize_records(frame: pd.DataFrame, columns) -> pd.DataFrame:
    """Return records in the canonical schema with whitespace-normalized values."""
    required = list(columns)
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError("The imported table is missing required columns: " + ", ".join(missing))
    normalized = frame.loc[:, required].copy()
    for column in required:
        normalized[column] = normalized[column].fillna("").astype(str).str.strip()
    return normalized


def record_id(row: pd.Series) -> str:
    identity = "\x1f".join(str(row.get(column, "")).strip().casefold() for column in RECORD_KEY_COLUMNS)
    return hashlib.sha256(identity.encode("utf-8")).hexdigest()[:16]


def with_record_ids(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result.insert(0, "record_id", result.apply(record_id, axis=1))
    return result


def merge_records(current: pd.DataFrame, incoming: pd.DataFrame, columns) -> tuple[pd.DataFrame, dict]:
    """Add new evidence and update matching records without duplicating them."""
    current = normalize_records(current, columns) if not current.empty else pd.DataFrame(columns=list(columns))
    incoming = normalize_records(incoming, columns)
    incoming = incoming.drop_duplicates(subset=list(RECORD_KEY_COLUMNS), keep="last")
    current_tagged = with_record_ids(current).drop_duplicates("record_id", keep="last")
    incoming_tagged = with_record_ids(incoming).drop_duplicates("record_id", keep="last")
    current_ids = current_tagged.set_index("record_id", drop=False) if not current.empty else None
    incoming_ids = incoming_tagged.set_index("record_id", drop=False)
    old_ids = set() if current_ids is None else set(current_ids.index)
    new_ids = set(incoming_ids.index)
    added = len(new_ids - old_ids)
    updated = unchanged = 0
    if current_ids is not None:
        for key in new_ids & old_ids:
            if current_ids.loc[key, list(columns)].tolist() == incoming_ids.loc[key, list(columns)].tolist():
                unchanged += 1
            else:
                updated += 1
    combined = pd.concat([current_tagged, incoming_tagged], ignore_index=True).drop_duplicates(
        subset=["record_id"], keep="last"
    )
    return combined.loc[:, list(columns)].reset_index(drop=True), {
        "added": added, "updated": updated, "unchanged": unchanged
    }


def delete_records(current: pd.DataFrame, record_ids) -> pd.DataFrame:
    ids = set(record_ids)
    if not ids or current.empty:
        return current.copy()
    tagged = with_record_ids(current)
    return tagged.loc[~tagged["record_id"].isin(ids), current.columns].reset_index(drop=True)
