#!/usr/bin/env python3
import os
import gzip
import tempfile
import shutil
import argparse
import time
import signal
import warnings
from typing import Optional, Tuple, Dict, List

import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.DSSP import DSSP


# Helix types + beta strand/bridge are "structured"
STRUCTURED_STATES = {"H", "G", "I", "E", "B"}


def read_ids(csv_path: str) -> List[str]:
    df = pd.read_csv(csv_path)
    if df.shape[1] == 0:
        raise ValueError(f"No columns found in {csv_path}")

    col = "pdb_id" if "pdb_id" in df.columns else df.columns[0]
    ids = (
        df[col]
        .astype(str)
        .str.strip()
        .str.lower()
        .replace({"nan": None})
        .dropna()
        .tolist()
    )
    return ids


def find_structure_file(structures_dir: str, pdb_id: str) -> Optional[str]:
    candidates = [
        os.path.join(structures_dir, f"{pdb_id}.cif.gz"),
        os.path.join(structures_dir, f"{pdb_id}.cif"),
        os.path.join(structures_dir, f"{pdb_id.upper()}.cif.gz"),
        os.path.join(structures_dir, f"{pdb_id.upper()}.cif"),
    ]
    for p in candidates:
        if os.path.exists(p):
            return p
    return None


def gunzip_to_temp(cif_gz_path: str) -> str:
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".cif")
    tmp_path = tmp.name
    tmp.close()
    with gzip.open(cif_gz_path, "rb") as f_in, open(tmp_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return tmp_path


def compute_max_internal_unstructured_run(ss_codes: List[str]) -> int:
    """
    Unstructured = DSSP code NOT in STRUCTURED_STATES.
    We ignore terminal runs (any run touching first or last residue).
    """
    n = len(ss_codes)
    if n == 0:
        return 0

    is_unstruct = [c not in STRUCTURED_STATES for c in ss_codes]

    max_run = 0
    i = 0
    while i < n:
        if not is_unstruct[i]:
            i += 1
            continue
        j = i
        while j < n and is_unstruct[j]:
            j += 1

        touches_left = (i == 0)
        touches_right = (j == n)
        if not touches_left and not touches_right:
            max_run = max(max_run, j - i)

        i = j

    return max_run


class _Timeout(Exception):
    pass


def _alarm_handler(signum, frame):
    raise _Timeout("DSSP timed out")


def run_dssp_on_cif(
    cif_path: str,
    pdb_id: str,
    dssp_exe: str,
    timeout_sec: int
) -> DSSP:
    """
    Runs DSSP via Biopython with a UNIX alarm timeout (works on macOS).
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, cif_path)
    model = next(structure.get_models())

    old_handler = signal.signal(signal.SIGALRM, _alarm_handler)
    try:
        signal.alarm(int(timeout_sec))
        dssp = DSSP(model, cif_path, dssp=dssp_exe)
        return dssp
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)


def analyze_entry(
    cif_or_gz_path: str,
    pdb_id: str,
    dssp_exe: str,
    ratio_thresh: float,
    timeout_sec: int,
    chain_select: str = "all"  # "all" or "top2"
) -> Tuple[bool, Dict]:
    """
    Returns (keep, details)
    Discard if ANY selected chain has max_internal_unstructured_run / chain_length > ratio_thresh
    """
    details = {"pdb_id": pdb_id, "status": None, "reason": None}
    tmp_cif = None

    try:
        cif_path = cif_or_gz_path
        if cif_or_gz_path.endswith(".gz"):
            tmp_cif = gunzip_to_temp(cif_or_gz_path)
            cif_path = tmp_cif

        # Optional: silence noisy mmCIF dictionary version warnings (doesn't affect most runs)
        warnings.filterwarnings(
            "ignore",
            message="The version in dictionary mmcif_pdbx.dic is lower than requested*",
            category=UserWarning,
        )

        dssp = run_dssp_on_cif(cif_path, pdb_id, dssp_exe, timeout_sec)

        # Collect DSSP codes per chain in residue order
        chain_map: Dict[str, List[Tuple[int, str, str]]] = {}
        for (chain_id, res_id) in dssp.keys():
            hetflag, resseq, icode = res_id
            code = dssp[(chain_id, res_id)][2]  # DSSP secondary structure code
            chain_map.setdefault(chain_id, []).append((resseq, str(icode), code))

        if not chain_map:
            details["status"] = "discarded"
            details["reason"] = "No DSSP-assigned residues found"
            return False, details

        # If a structure has more than 2 chains, you can choose to evaluate only the top2 longest
        chain_items = []
        for chain_id, triples in chain_map.items():
            triples.sort(key=lambda x: (x[0], x[1]))
            ss_codes = [t[2] for t in triples]
            chain_items.append((chain_id, ss_codes))

        if chain_select == "top2" and len(chain_items) > 2:
            chain_items.sort(key=lambda x: len(x[1]), reverse=True)
            chain_items = chain_items[:2]
            details["note"] = f"More than 2 chains present; evaluated top2 by DSSP length."

        per_chain_rows = []
        bad_chain = False

        for chain_id, ss_codes in chain_items:
            L = len(ss_codes)
            max_internal = compute_max_internal_unstructured_run(ss_codes)
            ratio = (max_internal / L) if L > 0 else 0.0
            per_chain_rows.append((chain_id, L, max_internal, ratio))
            if L > 0 and ratio > ratio_thresh:
                bad_chain = True

        worst = max(per_chain_rows, key=lambda x: x[3])
        details.update({
            "worst_chain": worst[0],
            "chain_length": worst[1],
            "max_internal_unstruct": worst[2],
            "max_internal_unstruct_ratio": worst[3],
            "all_chains": ";".join(
                f"{cid}:L={L},max={mx},ratio={rt:.3f}" for cid, L, mx, rt in per_chain_rows
            )
        })

        if bad_chain:
            details["status"] = "discarded"
            details["reason"] = f"Internal unstructured segment > {int(ratio_thresh*100)}% of chain length"
            return False, details

        details["status"] = "kept"
        details["reason"] = "Passed rigidity filter"
        return True, details

    except _Timeout:
        details["status"] = "discarded"
        details["reason"] = f"Error: DSSP timeout after {timeout_sec}s"
        return False, details

    except Exception as e:
        details["status"] = "discarded"
        details["reason"] = f"Error: {type(e).__name__}: {e}"
        return False, details

    finally:
        if tmp_cif and os.path.exists(tmp_cif):
            try:
                os.remove(tmp_cif)
            except OSError:
                pass


def write_outputs(base_dir: str, kept: List[str], discarded: List[str], report_rows: List[Dict],
                  out_kept: str, out_discarded: str, out_report: str) -> None:
    pd.DataFrame({"pdb_id": kept}).to_csv(os.path.join(base_dir, out_kept), index=False)
    pd.DataFrame({"pdb_id": discarded}).to_csv(os.path.join(base_dir, out_discarded), index=False)
    pd.DataFrame(report_rows).to_csv(os.path.join(base_dir, out_report), index=False)


def main():
    ap = argparse.ArgumentParser(description="Secondary-structure rigidity filter for heterodimers using DSSP.")
    ap.add_argument("--base_dir", required=True, help="Path to your postprocessing folder.")
    ap.add_argument("--ids_csv", default="primary_pdb_ids.csv", help="CSV containing pdb ids (default: primary_pdb_ids.csv).")
    ap.add_argument("--structures_dir", default="structures", help="Folder containing mmCIF files (default: structures).")
    ap.add_argument("--dssp_exe", default="mkdssp", help="DSSP executable name/path (default: mkdssp).")
    ap.add_argument("--ratio_thresh", type=float, default=0.30, help="Discard if internal unstructured run ratio exceeds this (default: 0.30).")
    ap.add_argument("--timeout_sec", type=int, default=120, help="Timeout (seconds) per entry for DSSP (default: 120).")
    ap.add_argument("--progress_every", type=int, default=25, help="Print progress every N entries (default: 25).")
    ap.add_argument("--slow_sec", type=int, default=30, help="Warn if a single entry takes > N seconds (default: 30).")
    ap.add_argument("--checkpoint_every", type=int, default=200, help="Write checkpoint outputs every N entries (default: 200).")
    ap.add_argument("--chain_select", choices=["all", "top2"], default="all",
                    help="If structure has >2 chains, evaluate 'all' or only 'top2' longest by DSSP length (default: all).")
    ap.add_argument("--out_kept", default="kept_pdb_ids.csv", help="Output CSV for kept ids.")
    ap.add_argument("--out_discarded", default="discarded_pdb_ids.csv", help="Output CSV for discarded ids.")
    ap.add_argument("--out_report", default="rigidity_filter_report.csv", help="Detailed report CSV.")
    args = ap.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    ids_csv_path = os.path.join(base_dir, args.ids_csv)
    structures_dir = os.path.join(base_dir, args.structures_dir)

    pdb_ids = read_ids(ids_csv_path)

    print(f"Base dir: {base_dir}", flush=True)
    print(f"IDs CSV: {ids_csv_path}", flush=True)
    print(f"Structures dir: {structures_dir}", flush=True)
    print(f"Total IDs read: {len(pdb_ids)}", flush=True)
    print(f"DSSP exe: {args.dssp_exe} | timeout: {args.timeout_sec}s", flush=True)

    kept: List[str] = []
    discarded: List[str] = []
    report_rows: List[Dict] = []

    t0 = time.time()

    try:
        for idx, pdb_id in enumerate(pdb_ids, start=1):
            entry_t0 = time.time()

            if idx == 1 or idx % args.progress_every == 0:
                elapsed = time.time() - t0
                rate = idx / elapsed if elapsed > 0 else 0.0
                print(
                    f"[{idx}/{len(pdb_ids)}] starting {pdb_id} "
                    f"(kept={len(kept)} discarded={len(discarded)} rate={rate:.2f}/s)",
                    flush=True
                )

            path = find_structure_file(structures_dir, pdb_id)
            if path is None:
                details = {
                    "pdb_id": pdb_id,
                    "status": "discarded",
                    "reason": "Structure file not found",
                    "worst_chain": None,
                    "chain_length": None,
                    "max_internal_unstruct": None,
                    "max_internal_unstruct_ratio": None,
                    "all_chains": None
                }
                discarded.append(pdb_id)
                report_rows.append(details)
                continue

            keep, details = analyze_entry(
                path, pdb_id,
                dssp_exe=args.dssp_exe,
                ratio_thresh=args.ratio_thresh,
                timeout_sec=args.timeout_sec,
                chain_select=args.chain_select
            )

            report_rows.append(details)
            if keep:
                kept.append(pdb_id)
            else:
                discarded.append(pdb_id)

            entry_elapsed = time.time() - entry_t0
            if entry_elapsed > args.slow_sec:
                print(f"  ⚠️ slow: {pdb_id} took {entry_elapsed:.1f}s ({details.get('reason')})", flush=True)

            # Checkpoint writing so you never lose progress
            if idx % args.checkpoint_every == 0:
                write_outputs(base_dir, kept, discarded, report_rows,
                              args.out_kept, args.out_discarded, args.out_report)
                print(f"  💾 checkpoint saved at {idx}", flush=True)

    except KeyboardInterrupt:
        print("\nKeyboardInterrupt: saving checkpoint before exit...", flush=True)
        write_outputs(base_dir, kept, discarded, report_rows,
                      args.out_kept, args.out_discarded, args.out_report)
        print("Checkpoint saved. Exiting.", flush=True)
        return

    # Final write
    write_outputs(base_dir, kept, discarded, report_rows,
                  args.out_kept, args.out_discarded, args.out_report)

    print("Done.", flush=True)
    print(f"Kept: {len(kept)} -> {os.path.join(base_dir, args.out_kept)}", flush=True)
    print(f"Discarded: {len(discarded)} -> {os.path.join(base_dir, args.out_discarded)}", flush=True)
    print(f"Report: {os.path.join(base_dir, args.out_report)}", flush=True)


if __name__ == "__main__":
    main()
