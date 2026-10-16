"""Read S7 Expected columns only, after sealing all eleven predictions.

The literal table remains primary. A separate, systematic frozen-taxonomy view
is diagnostic; it never modifies predictions, model parameters or thresholds.
"""
import argparse
import csv
from datetime import datetime, timezone
import json
import re
from pathlib import Path

import openpyxl

from .export_real_predictions import sha
from .real_truth import verify_complete_batch, xml_sheet, taxonomy_index, resolve_taxonomy


TREE_SHA = "ff2cdf88a6f46da1f4155b36a4d9a96a838d415f9d609b8c8a284019f6e11b20"


def expected_rows(cells):
    for cell, value in (("B3", "Library ID"), ("D3", "Species"),
                        ("I3", "Expected Orders"), ("J3", "Expected Families")):
        if cells.get(cell) != value:
            raise ValueError("Unexpected S7 truth layout: " + cell)
    result = []
    row_numbers = sorted({int(re.search(r"\d+$", key).group()) for key in cells})
    for row in row_numbers:
        if row < 4 or not cells.get(f"B{row}"):
            continue
        library = str(cells[f"B{row}"])
        match = re.fullmatch(r"EM-(\d+)", library)
        if not match or not 1 <= int(match[1]) <= 11:
            raise ValueError("Unexpected library ID: " + library)
        record = {"sample": f"EMseq{int(match[1]):03d}", "library_id": library,
            "source_cell": f"D{row}", "library_source_cell": f"B{row}",
            "species": str(cells[f"D{row}"]).strip(), "family": str(cells[f"J{row}"]).strip(),
            "order": str(cells[f"I{row}"]).strip(), "family_source_cell": f"J{row}",
            "order_source_cell": f"I{row}", "mapping_rule": "explicit_S7_Expected_annotation",
            "mix_name": str(cells.get(f"A{row}", "")), "component_number": int(cells[f"C{row}"]),
            "weight_g": float(cells[f"E{row}"]), "tissue_id": str(cells.get(f"F{row}", "")),
            "tissue_type": str(cells.get(f"G{row}", "")), "proportion_of_mix_pct": float(cells[f"H{row}"])}
        if record["weight_g"] <= 0 or record["family"] in ("", ".", "None") or record["order"] in ("", ".", "None"):
            raise ValueError("Incomplete or nonpositive S7 component")
        result.append(record)
    if {r["sample"] for r in result} != {f"EMseq{i:03d}" for i in range(1, 12)}:
        raise ValueError("S7 must cover all eleven libraries")
    return result


def write_rows(path, rows):
    with Path(path).open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse(batch, workbook, taxonomy, output):
    batch, workbook, taxonomy, output = map(Path, (batch, workbook, taxonomy, output))
    verify_complete_batch(batch)
    seal = json.loads((output / "PREDICTION_EXPORT.json").read_text())
    for name, digest in seal["artifact_sha256"].items():
        if sha(output / name) != digest:
            raise ValueError("Sealed prediction changed: " + name)
    protocol = json.loads((batch / "PROTOCOL.json").read_text())
    if sha(workbook) != protocol["truth_workbook_sha256"] or sha(taxonomy) != TREE_SHA:
        raise ValueError("Truth workbook or frozen taxonomy checksum mismatch")
    started = datetime.now(timezone.utc).isoformat()
    wb = openpyxl.load_workbook(workbook, read_only=True, data_only=True)
    ws = wb["Table S7"]
    cells = {c.coordinate: c.value for row in ws for c in row if c.value is not None}
    xml = xml_sheet(workbook, "Table S7")
    for key in set(cells) | set(xml):
        value, raw = cells.get(key), xml.get(key, "")
        same = value == float(raw) if isinstance(value, (int, float)) else str(value or "") == raw
        if not same:
            raise ValueError("openpyxl/XML mismatch: " + key)
    primary, independent = expected_rows(cells), expected_rows(xml)
    keys = ("sample", "species", "family", "order", "source_cell", "weight_g", "proportion_of_mix_pct")
    if [tuple(r[k] for k in keys) for r in primary] != [tuple(r[k] for k in keys) for r in independent]:
        raise ValueError("Independent S7 component extraction mismatch")
    species, genera = taxonomy_index(taxonomy)
    canonical, disagreements = [], []
    for row in primary:
        try:
            mapped = resolve_taxonomy(row["species"], species, genera)
        except ValueError as error:
            # Keep the explicit annotation, never fill it from model predictions.
            mapped = {"family": row["family"], "order": row["order"],
                      "mapping_rule": "S7_annotation_retained_unresolved_frozen_taxonomy: " + str(error)}
        canonical.append({**row, **mapped})
        if any(row[k] != mapped[k] for k in ("family", "order")):
            disagreements.append({"sample": row["sample"], "source_cell": row["source_cell"], "species": row["species"],
                "S7_family": row["family"], "frozen_Kew_family": mapped["family"],
                "S7_order": row["order"], "frozen_Kew_order": mapped["order"], "mapping_rule": mapped["mapping_rule"]})
    write_rows(output / "truth_mapping.csv", primary)
    write_rows(output / "truth_mapping_kew.csv", canonical)
    if disagreements:
        write_rows(output / "taxonomy_disagreements.csv", disagreements)
    audit_path = output / "S7_SOURCE_AUDIT.json"
    audit = json.loads(audit_path.read_text()) if audit_path.exists() else {}
    audit.update({"parser_audit_started_utc": started, "prediction_sealed_utc": seal["sealed_utc"],
        "sealed_before_truth": seal["sealed_utc"] < started, "workbook_sha256": sha(workbook),
        "taxonomy_sha256": sha(taxonomy), "sheet": "Table S7", "source_range": "B4:J67",
        "primary_columns": "B=Library ID, D=Species, I=Expected Orders, J=Expected Families",
        "legacy_predicted_columns_K_L_used_as_truth": False, "all_text_cells_match_exactly": True,
        "all_numeric_cells_match_as_binary64": True, "nonempty_cells_checked": len(cells),
        "component_rows_independently_checked": len(primary), "libraries": 11,
        "literal_truth_kept_unchanged": True, "taxonomy_disagreements": disagreements,
        "secondary_mapping_applied_to_every_component": True,
        "secondary_unresolved_rows": sum(r["mapping_rule"].startswith("S7_annotation_retained") for r in canonical)})
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    wb.close()
    if sha(workbook) != protocol["truth_workbook_sha256"]:
        raise ValueError("Workbook changed during read-only extraction")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("batch", "workbook", "taxonomy", "output"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    parse(args.batch, args.workbook, args.taxonomy, args.output)
