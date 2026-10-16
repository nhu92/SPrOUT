"""Independent ZIP/XML truth read and direct candidate-list count audit."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile


def xml_truth(path):
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = ["".join(t.text or "" for t in item.findall(".//s:t", ns)) for item in root]
        book = ET.fromstring(archive.read("xl/workbook.xml"))
        sheet = next(s for s in book.findall("s:sheets/s:sheet", ns) if s.attrib["name"] == "Table S1")
        rid = sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        target = next(r.attrib["Target"] for r in relationships if r.attrib["Id"] == rid).lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        root = ET.fromstring(archive.read(target))
        truth = {}
        for row in root.findall("s:sheetData/s:row", ns):
            if int(row.attrib["r"]) < 4:
                continue
            values = {}
            for cell in row.findall("s:c", ns):
                column = "".join(c for c in cell.attrib["r"] if c.isalpha())
                if column not in {"A", "H"}:
                    continue
                value = cell.find("s:v", ns)
                if cell.attrib.get("t") == "inlineStr":
                    text = "".join(t.text or "" for t in cell.findall(".//s:t", ns))
                elif value is None:
                    continue
                elif cell.attrib.get("t") == "s":
                    text = shared[int(value.text)]
                else:
                    text = value.text
                values[column] = text
            try:
                identity = int(float(values.get("A", "nan")))
            except (ValueError, OverflowError):
                continue
            if 1 <= identity <= 30 and values.get("H"):
                truth[identity] = values["H"].strip()
        assert set(truth) == set(range(1, 31))
        return truth


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def audit(batch, evaluation, truth_file):
    batch, evaluation = Path(batch), Path(evaluation)
    assert (batch / "PREDICTIONS_COMPLETE").exists()
    truth = xml_truth(truth_file)
    samples = json.loads((batch / "selection.json").read_text())["selected"]
    rows = read_csv(evaluation / "by_sample.csv")
    detail = read_csv(evaluation / "by_family.csv")
    universe = {r["family"] for r in detail}
    recomputed = []
    for row in rows:
        sample = next(s for s in samples if s["sample"] == row["sample"])
        correction = {1: "Cabombaceae", 25: "Capparaceae"} if row["view"] == "harmonized_two_taxa" else {}
        expected = {correction.get(i, truth[i]) for i in map(int, row["sample"].split("x"))}
        called = set((batch / "samples" / sample["project"] / (row["method"] + ".txt")).read_text().splitlines())
        counts = dict(TP=len(expected & called), FP=len(called - expected), FN=len(expected - called),
                      TN=len(universe - expected - called))
        assert all(counts[k] == int(row[k]) for k in counts), (row["sample"], row["method"], row["view"])
        assert expected == set(row["expected"].split(";"))
        assert called == set(filter(None, row["predicted"].split(";")))
        recomputed.append({**row, **counts})
    summaries = read_csv(evaluation / "summary.csv")
    for total in summaries:
        group = total["group"]
        sub = [r for r in recomputed if r["view"] == total["view"] and r["method"] == total["method"] and
            (group == "all100" or group == "unique98" and r["first_composition"] == "True" or
             group == "historical30" and r["historical30"] == "True" or group == "remaining70" and r["historical30"] == "False")]
        assert len(sub) == int(total["samples"])
        assert all(sum(r[k] for r in sub) == int(total[k]) for k in ["TP", "FP", "FN", "TN"])
    result = {"independent_truth_parser": "stdlib ZIP/XML, cells A and H of Table S1",
              "truth_sha256": hashlib.sha256(Path(truth_file).read_bytes()).hexdigest(),
              "samples": len(samples), "sample_method_view_rows_checked": len(rows),
              "summary_rows_checked": len(summaries), "fixed_family_universe": len(universe),
              "all_candidate_lists_and_counts_match": True}
    (evaluation / "COUNT_AUDIT.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", required=True)
    parser.add_argument("--evaluation", required=True)
    parser.add_argument("--truth", required=True)
    args = parser.parse_args()
    audit(args.batch, args.evaluation, args.truth)
