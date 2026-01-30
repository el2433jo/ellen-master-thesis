#!/usr/bin/env python3
"""
Split Complex Portal JSON into two CSV files:

1) heterodimers where BOTH participants are proteins and each chain length <= MAX_LEN
2) everything else

Each CSV row includes FASTA-formatted sequences (one FASTA record per chain).

Progress report is printed at the end (counts + output paths).

Default paths:
  in:  /Users/ellenjonsson/Documents/Ellen's scripts/postprocessing/yeast/yeast_complexportal.json
  out: /Users/ellenjonsson/Documents/Ellen's scripts/postprocessing/yeast/

Outputs:
  - yeast_heterodimers_le{MAX_LEN}.csv
  - yeast_other_complexes.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

UNIPROTKB_ID_RE = re.compile(r"^uniprotkb_[A-Za-z0-9]+(?:-\d+)?$")  # allow isoform-like suffixes


def safe_get(d: Dict[str, Any], path: Tuple[str, ...]) -> Any:
    cur: Any = d
    for p in path:
        if not isinstance(cur, dict) or p not in cur:
            return None
        cur = cur[p]
    return cur


def to_uniprot(interactor_id: str) -> str:
    return interactor_id.split("_", 1)[1] if interactor_id.startswith("uniprotkb_") else interactor_id


def wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def make_fasta(header: str, seq: Optional[str]) -> str:
    if not seq:
        return ""
    return f">{header}\n{wrap_fasta(seq)}"


def build_protein_index(data: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """
    Returns mapping interactor_id -> dict with uniprot,label,seq,len,organism info.
    """
    idx: Dict[str, Dict[str, Any]] = {}
    for obj in data:
        if not isinstance(obj, dict) or obj.get("object") != "interactor":
            continue

        interactor_id = obj.get("id")
        if not isinstance(interactor_id, str):
            continue

        # Identify proteins
        type_name = safe_get(obj, ("type", "name"))
        type_id = safe_get(obj, ("type", "id"))
        is_protein = (type_name == "protein") or (type_id == "MI:0326")
        if not is_protein:
            continue

        seq = obj.get("sequence") if isinstance(obj.get("sequence"), str) else None
        idx[interactor_id] = {
            "interactor_id": interactor_id,
            "uniprot": to_uniprot(interactor_id),
            "label": obj.get("label") if isinstance(obj.get("label"), str) else "",
            "sequence": seq or "",
            "length": len(seq) if seq else 0,
            "taxid": str(safe_get(obj, ("organism", "taxid")) or ""),
            "scientific": str(safe_get(obj, ("organism", "scientific")) or ""),
        }
    return idx


def extract_from_participant_like(p: Dict[str, Any]) -> Tuple[List[str], Optional[int]]:
    ids: List[str] = []

    for key in ("interactorRef", "interactor_id", "interactorId"):
        v = p.get(key)
        if isinstance(v, str) and UNIPROTKB_ID_RE.match(v):
            ids.append(v)

    interactor_obj = p.get("interactor")
    if isinstance(interactor_obj, dict):
        vid = interactor_obj.get("id")
        if isinstance(vid, str) and UNIPROTKB_ID_RE.match(vid):
            ids.append(vid)

    identifiers = p.get("identifiers")
    if isinstance(identifiers, list):
        for item in identifiers:
            if isinstance(item, dict):
                if item.get("db") == "uniprotkb" and isinstance(item.get("id"), str):
                    ids.append(f"uniprotkb_{item['id']}")

    stoich: Optional[int] = None
    for sk in ("stoichiometry", "stoichiometryValue"):
        sv = p.get(sk)
        if isinstance(sv, int):
            stoich = sv
        elif isinstance(sv, str) and sv.isdigit():
            stoich = int(sv)

    return ids, stoich


def recursive_collect_uniprotkb_ids(x: Any) -> List[str]:
    found: List[str] = []
    if isinstance(x, dict):
        for v in x.values():
            found.extend(recursive_collect_uniprotkb_ids(v))
    elif isinstance(x, list):
        for v in x:
            found.extend(recursive_collect_uniprotkb_ids(v))
    elif isinstance(x, str):
        if UNIPROTKB_ID_RE.match(x):
            found.append(x)
    return found


def extract_participants(interaction: Dict[str, Any]) -> Tuple[List[str], Optional[Counter]]:
    participant_lists: List[List[Any]] = []
    for key in ("participants", "interactionParticipants", "interactionComponents", "components"):
        v = interaction.get(key)
        if isinstance(v, list):
            participant_lists.append(v)

    ids: List[str] = []
    stoich_counter: Optional[Counter] = None

    if participant_lists:
        counts = Counter()
        any_stoich = False

        for plist in participant_lists:
            for item in plist:
                if not isinstance(item, dict):
                    continue
                item_ids, stoich = extract_from_participant_like(item)
                if not item_ids:
                    continue
                for iid in item_ids:
                    if stoich is not None:
                        any_stoich = True
                        counts[iid] += stoich
                    else:
                        counts[iid] += 1
                    ids.append(iid)

        if counts:
            stoich_counter = counts if any_stoich else Counter(ids)

        uniq = []
        seen = set()
        for iid in ids:
            if iid not in seen:
                uniq.append(iid)
                seen.add(iid)
        return uniq, stoich_counter

    rec = recursive_collect_uniprotkb_ids(interaction)
    if rec:
        return sorted(set(rec)), Counter(rec)

    return [], None


def classify(
    interaction: Dict[str, Any],
    protein_index: Dict[str, Dict[str, Any]],
    max_len: int,
) -> Tuple[bool, List[str], List[Dict[str, Any]]]:
    """
    Returns (is_heterodimer_lemax, reasons_if_not, participants_protein_dicts)
    """
    participant_ids, counter = extract_participants(interaction)
    reasons: List[str] = []

    if not participant_ids:
        reasons.append("No participant UniProtKB IDs found.")

    participants: List[Dict[str, Any]] = []
    missing: List[str] = []
    for pid in participant_ids:
        p = protein_index.get(pid)
        if p is None:
            missing.append(pid)
        else:
            participants.append(p)

    if missing:
        reasons.append("Missing protein records for: " + ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else ""))

    uniq_uniprot = sorted({p["uniprot"] for p in participants})
    if len(uniq_uniprot) != 2:
        reasons.append(f"Not a heterodimer (need 2 distinct protein participants; found {len(uniq_uniprot)}).")

    if counter is not None and len(counter) > 0 and len(uniq_uniprot) == 2:
        # If we can, enforce 1:1 total
        prot_ids = [p["interactor_id"] for p in participants]
        total = sum(counter.get(iid, 0) for iid in prot_ids)
        if total != 2:
            reasons.append(f"Counts/stoichiometry not consistent with 1:1 dimer (total={total}).")

    for p in participants:
        if not p["sequence"]:
            reasons.append(f"Missing sequence for {p['uniprot']}.")
        elif p["length"] > max_len:
            reasons.append(f"Chain {p['uniprot']} length {p['length']} > {max_len}.")

    ok = (
        len(uniq_uniprot) == 2
        and not any("Missing sequence" in r for r in reasons)
        and not any("length" in r and ">" in r for r in reasons)
        and not any(r.startswith("Not a heterodimer") for r in reasons)
        and not any("Counts/stoichiometry" in r for r in reasons)
    )

    return ok, reasons, participants


def main() -> None:
    default_dir = Path(r"/Users/ellenjonsson/Documents/Ellen's scripts/postprocessing/yeast")

    ap = argparse.ArgumentParser(description="Split Complex Portal JSON into CSVs with FASTA sequences.")
    ap.add_argument("--input", type=Path, default=default_dir / "yeast_complexportal.json",
                    help="Path to Complex Portal JSON file.")
    ap.add_argument("--outdir", type=Path, default=default_dir, help="Output directory.")
    ap.add_argument("--max-len", type=int, default=500, help="Max amino acids per chain.")
    args = ap.parse_args()

    inpath: Path = args.input
    outdir: Path = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    with inpath.open("r", encoding="utf-8") as f:
        root = json.load(f)

    data = root.get("data")
    if not isinstance(data, list):
        raise ValueError("Input JSON does not contain a top-level 'data' list.")

    protein_index = build_protein_index(data)

    interactions: List[Dict[str, Any]] = [
        obj for obj in data
        if isinstance(obj, dict) and obj.get("object") == "interaction"
    ]

    hetero_rows: List[Dict[str, Any]] = []
    other_rows: List[Dict[str, Any]] = []

    for inter in interactions:
        complex_id = inter.get("id") if isinstance(inter.get("id"), str) else "(no id)"
        ok, reasons, participants = classify(inter, protein_index, max_len=args.max_len)

        # Try to pick the two chains for output (best-effort)
        # Sort by uniprot to be stable
        participants_sorted = sorted(participants, key=lambda p: p["uniprot"])
        chainA = participants_sorted[0] if len(participants_sorted) >= 1 else None
        chainB = participants_sorted[1] if len(participants_sorted) >= 2 else None

        def chain_fields(prefix: str, p: Optional[Dict[str, Any]]) -> Dict[str, Any]:
            if not p:
                return {
                    f"{prefix}_uniprot": "",
                    f"{prefix}_label": "",
                    f"{prefix}_length": "",
                    f"{prefix}_fasta": "",
                }
            header = f"{p['uniprot']}|{p.get('label','')}".strip("|")
            return {
                f"{prefix}_uniprot": p["uniprot"],
                f"{prefix}_label": p.get("label", ""),
                f"{prefix}_length": p.get("length", ""),
                f"{prefix}_fasta": make_fasta(header, p.get("sequence", "")),
            }

        row = {
            "complex_id": complex_id,
            "organism_taxid": str(safe_get(inter, ("organism", "taxid")) or ""),
            "organism_scientific": str(safe_get(inter, ("organism", "scientific")) or ""),
            "participant_count_extracted": len(participants),
            "passes_heterodimer_le_maxlen": "yes" if ok else "no",
            "failure_reasons": "; ".join(reasons) if reasons else "",
        }
        row.update(chain_fields("chainA", chainA))
        row.update(chain_fields("chainB", chainB))

        if ok:
            hetero_rows.append(row)
        else:
            other_rows.append(row)

    hetero_path = outdir / f"yeast_heterodimers_le{args.max_len}.csv"
    other_path = outdir / "yeast_other_complexes.csv"

    fieldnames = [
        "complex_id",
        "organism_taxid",
        "organism_scientific",
        "participant_count_extracted",
        "passes_heterodimer_le_maxlen",
        "failure_reasons",
        "chainA_uniprot", "chainA_label", "chainA_length", "chainA_fasta",
        "chainB_uniprot", "chainB_label", "chainB_length", "chainB_fasta",
    ]

    def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
        with path.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r in rows:
                w.writerow(r)

    write_csv(hetero_path, hetero_rows)
    write_csv(other_path, other_rows)

    # Progress report / summary (printed after the run)
    print("\n=== Complex Portal split summary ===")
    print(f"Input file: {inpath}")
    print(f"Output dir: {outdir}")
    print(f"Protein interactors indexed: {len(protein_index)}")
    print(f"Total interaction objects:   {len(interactions)}")
    print(f"Heterodimers <= {args.max_len} aa: {len(hetero_rows)}")
    print(f"Other complexes:             {len(other_rows)}")
    print(f"Wrote: {hetero_path}")
    print(f"Wrote: {other_path}")
    print("===================================\n")


if __name__ == "__main__":
    main()