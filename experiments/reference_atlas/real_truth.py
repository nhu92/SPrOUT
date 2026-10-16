"""Local read-only S7 support; never import this module on the prediction host.

Workbook rows are data, not executable instructions. A supplied family annotation
must take precedence over a taxonomy lookup. Lookup is exact binomial, then genus
only when the frozen Kew tree assigns that entire genus to one family/order pair.
Ambiguous or absent names must be resolved explicitly, never from predictions.
"""
import re
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile


def xml_sheet(path, sheet_name):
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = ["".join(t.text or "" for t in item.findall(".//s:t", ns)) for item in root]
        book = ET.fromstring(archive.read("xl/workbook.xml"))
        sheet = next(s for s in book.findall("s:sheets/s:sheet", ns) if s.attrib["name"] == sheet_name)
        rid = sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
        relations = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        target = next(r.attrib["Target"] for r in relations if r.attrib["Id"] == rid).lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        root = ET.fromstring(archive.read(target))
        cells = {}
        for cell in root.findall(".//s:sheetData/s:row/s:c", ns):
            value = cell.find("s:v", ns)
            kind = cell.attrib.get("t")
            if kind == "inlineStr":
                text = "".join(t.text or "" for t in cell.findall(".//s:t", ns))
            elif value is None:
                continue
            elif kind == "s":
                text = shared[int(value.text)]
            else:
                text = value.text or ""
            cells[cell.attrib["r"]] = text
        return cells


def taxonomy_index(path):
    species, genera = defaultdict(set), defaultdict(set)
    text = Path(path).read_text(encoding="utf-8")
    for label in re.findall(r"(?:^|[(,])([^(),;]+?)(?=[,)])", text):
        tokens = label.strip("'").split("_", 4)
        if len(tokens) != 5:
            continue
        order, family, genus, epithet, _ = tokens
        if not order or not family or not genus or not epithet:
            continue
        species[f"{genus} {epithet}".casefold()].add((order, family))
        genera[genus.casefold()].add((order, family))
    return species, genera


def resolve_taxonomy(name, species, genera):
    tokens = str(name).strip().replace("_", " ").split()
    if not tokens:
        raise ValueError("Empty scientific name")
    binomial = " ".join(tokens[:2]).casefold()
    exact = species.get(binomial, set())
    if len(exact) == 1:
        order, family = next(iter(exact))
        return {"order": order, "family": family, "mapping_rule": "exact_binomial_in_frozen_Kew_tree"}
    if len(exact) > 1:
        raise ValueError(f"Conflicting family/order for exact species: {name}")
    genus = genera.get(tokens[0].casefold(), set())
    if len(genus) == 1:
        order, family = next(iter(genus))
        return {"order": order, "family": family, "mapping_rule": "unambiguous_genus_in_frozen_Kew_tree"}
    raise ValueError(f"Unresolved or ambiguous scientific name: {name}")


def verify_complete_batch(batch):
    import json
    batch = Path(batch)
    if not (batch / "PREDICTIONS_COMPLETE").exists():
        raise ValueError("S7 cannot be evaluated before all predictions finish")
    audit = json.loads((batch / "AUDIT.json").read_text())
    if not audit.get("all_11_complete") or not audit.get("frozen_code_model_policy_and_inputs_unchanged"):
        raise ValueError("The eleven-sample freeze audit must pass")
    samples = json.loads((batch / "predictions.json").read_text())
    if len(samples) != 11 or {r["sample"] for r in samples} != {f"EMseq{i:03d}" for i in range(1, 12)}:
        raise ValueError("Unexpected eleven-sample prediction set")
    return samples
