"""
Pipeline overview (what I do in this script)

1) I combine multiple RCSB custom report CSV exports (e.g. rcsb_report_1.csv ... rcsb_report_5.csv)
   into one clean CSV with a single header row.

2) I parse the combined report and reconstruct, per PDB ID:
   - the two protein entity sequences (Entity ID 1 and 2)
   - quality metadata (EM resolution or X-ray resolution limit)
   - ligand annotations (distinct ligand names and/or ligand counts)

   Important: RCSB's export often writes "continuation rows" where the PDB ID is blank.
   I handle that by carrying forward the last seen PDB ID until a new one appears.

3) I cluster by *exact heterodimer sequence pair* (order-independent):
      key = {sequence A, sequence B}
   This removes duplicates of the same complex solved under different conditions.

4) I choose ONE representative per exact-sequence cluster.
   By default I prefer:
      (a) fewer ligands
      (b) better resolution (lower is better)
   You can change that rule easily (see `representative_score`).

5) I remove "point-mutant duplicates" among representatives only:
   if two representatives differ by exactly 1 amino-acid substitution in one partner
   (and the other partner is identical), I keep only the better one (same scoring rule).

6) I write three outputs:
   - primary_pdb_ids.csv                (just one column: PDB ID)
   - exact_sequence_clusters.csv        (rep + all members, with ligand info)
   - point_mutant_decisions.csv         (what got removed in the point-mutant step)
"""

import csv
import glob
import hashlib
from collections import defaultdict
from itertools import combinations
from math import inf

# =========================
# 0) USER SETTINGS
# =========================

# Put your chunked downloads in the same folder and name them like:
#   rcsb_report_1.csv, rcsb_report_2.csv, ...
INPUT_GLOB = "rcsb_report_*.csv"

# Combined file I generate (single header row, all rows appended)
COMBINED_CSV = "rcsb_report_all.csv"

# Outputs
OUT_PRIMARY = "primary_pdb_ids.csv"
OUT_CLUSTERS = "exact_sequence_clusters.csv"
OUT_POINT_MUT = "point_mutant_decisions.csv"


# =========================
# 1) SMALL HELPER FUNCTIONS
# =========================

def to_float(x):
    """Convert a string to float, returning None if blank/unparseable."""
    try:
        return float(x) if x not in ("", None) else None
    except ValueError:
        return None

def hamming(a, b):
    """
    Hamming distance for equal-length strings.
    Returns None if lengths differ (so we don't treat indels as "point mutants").
    """
    if len(a) != len(b):
        return None
    return sum(x != y for x, y in zip(a, b))

def dimer_key(seq1, seq2):
    """
    Order-independent key for a heterodimer:
      A-B is the same as B-A
    """
    return tuple(sorted((seq1, seq2)))

def cluster_id_from_key(key):
    """
    Short stable ID for a cluster so it's readable in CSV.
    I hash the two sequences (order-independent by construction).
    """
    h = hashlib.sha1((key[0] + "|" + key[1]).encode("utf-8")).hexdigest()
    return h[:12]

def quality_score(meta):
    """
    Lower is better.
    If EM resolution exists, I use it.
    Otherwise I use X-ray high resolution limit.
    If neither exists, I treat it as 'worst' (infinity).
    """
    if meta.get("em_res") is not None:
        return meta["em_res"]
    if meta.get("xray_res") is not None:
        return meta["xray_res"]
    return inf

def representative_score(meta):
    """
    This is the rule I use to pick the 'best' representative.
    Lower tuple is better.

    Current policy:
      1) fewer ligands is better
      2) better resolution is better (lower Å)
    """
    ligand_count = meta.get("ligand_count")
    if ligand_count is None:
        ligand_count = inf  # if missing, treat as unknown/worst

    return (ligand_count, quality_score(meta))


# =========================
# 2) COMBINE MULTIPLE CSV CHUNKS INTO ONE
# =========================

def detect_real_header_and_rows(path):
    """
    RCSB exports often look like:
      line 1: "Identifier,StructureData,..."
      line 2: REAL HEADER ("Entry ID, EM Resolution..., Sequence, Entity ID, Ligand Name, ...")
      line 3+: data rows, including continuation rows (blank PDB ID)

    This function returns:
      header: the REAL header list
      rows: iterator over the data rows (lists)
    """
    f = open(path, "r", newline="")
    reader = csv.reader(f)

    line1 = next(reader, None)
    if line1 is None:
        f.close()
        return None, None

    line2 = next(reader, None)

    # If first line is the group header "Identifier,...", the second line is the real header
    if line1 and line1[0].strip() == "Identifier" and line2:
        header = line2
    else:
        # Otherwise we assume line1 itself is the real header
        header = line1
        # and line2 is actually the first data row
        # We'll handle that by yielding line2 as a data row below
        # (as long as it exists).
        # We'll keep reader positioned after line2.
        pass

    # Clean trailing empty header cells (RCSB sometimes leaves commas at the end)
    while header and header[-1] == "":
        header.pop()

    def row_iter():
        # If line1 was real header, then line2 is first data row (if present)
        if not (line1 and line1[0].strip() == "Identifier") and line2:
            yield line2
        for r in reader:
            yield r

    return header, row_iter(), f  # return f so caller can close it


def combine_reports(input_glob, combined_path):
    """
    Combine all chunked CSV files into a single CSV with a single header row.
    """
    paths = sorted(glob.glob(input_glob))
    if not paths:
        raise FileNotFoundError(
            f"No input files matched {input_glob}. "
            f"Make sure you named them like rcsb_report_1.csv, rcsb_report_2.csv, ..."
        )

    combined_header = None

    with open(combined_path, "w", newline="") as out_f:
        writer = None

        for i, path in enumerate(paths, start=1):
            header, rows, handle = detect_real_header_and_rows(path)
            if header is None:
                continue

            try:
                # On the first file, I set the header that will be used for all files
                if combined_header is None:
                    combined_header = header
                    writer = csv.writer(out_f)
                    writer.writerow(combined_header)
                else:
                    # Sanity check: header consistency across chunks
                    if header != combined_header:
                        raise ValueError(
                            f"Header mismatch in {path}.\n"
                            f"Expected: {combined_header}\n"
                            f"Found:    {header}\n"
                            "Make sure all chunks were generated with the same custom report columns."
                        )

                # Write all data rows from this chunk
                for r in rows:
                    if not r:
                        continue
                    writer.writerow(r)

            finally:
                handle.close()

    return paths


# =========================
# 3) PARSE THE COMBINED CSV INTO PER-PDB RECORDS
# =========================

def parse_combined_report(path):
    """
    Reads the combined CSV and reconstructs per PDB ID:
      - sequences for entity 1 and 2
      - method + resolution metadata
      - ligand info (ligand names and "Number of Distinct Non-polymer Entities" when provided)

    Output format:
      per_pdb[pdb_id] = {
         "seqs": { entity_id_str: sequence_str, ... },
         "meta": {...},
         "ligands": set([...])
      }
    """
    per_pdb = defaultdict(lambda: {"seqs": {}, "meta": {}, "ligands": set()})

    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)

        last_pdb = None
        last_meta = {"method": "", "em_res": None, "xray_res": None, "ligand_count": None}

        for rec in reader:
            # Primary identifiers
            pdb = (rec.get("PDB ID") or rec.get("Entry ID") or "").strip()

            # Quality metadata
            method = (rec.get("Experimental Method") or "").strip()
            em_res = to_float(rec.get("EM Resolution (Å)"))
            xray_res = to_float(rec.get("High Resolution Limit"))

            # Ligand metadata (these fields exist in your new file)
            ligand_count = to_float(rec.get("Number of Distinct Non-polymer Entities"))
            # ligand_count is conceptually an integer, but float parsing is fine (e.g. "3" -> 3.0)
            if ligand_count is not None:
                ligand_count = int(ligand_count)

            ligand_name = (rec.get("Ligand Name") or "").strip()

            # Protein entity data
            seq = (rec.get("Sequence") or "").strip()
            ent = (rec.get("Entity ID") or "").strip()

            # Handle continuation rows:
            # Many rows omit PDB ID and metadata, but still belong to the previous PDB entry.
            if pdb:
                last_pdb = pdb
                last_meta = {
                    "method": method,
                    "em_res": em_res,
                    "xray_res": xray_res,
                    "ligand_count": ligand_count,
                }
            else:
                pdb = last_pdb  # carry forward

            if not pdb:
                # If the file starts with blank rows for some reason, I just skip them.
                continue

            # Store metadata (safe to overwrite — it's the same per PDB ID)
            per_pdb[pdb]["meta"] = last_meta

            # Store ligand names (can appear on multiple rows)
            if ligand_name:
                per_pdb[pdb]["ligands"].add(ligand_name)

            # Store entity sequences (I expect exactly two entities for your filtered dataset)
            if seq and ent:
                per_pdb[pdb]["seqs"][ent] = seq

    return per_pdb


# =========================
# 4) CLUSTER BY EXACT SEQUENCE PAIR + CHOOSE REPRESENTATIVES
# =========================

def build_dimers(per_pdb):
    """
    Convert per_pdb into a simpler dict:
      dimers[pdb] = (seqA, seqB, meta, ligands_set)

    I skip entries that don't have exactly 2 sequences,
    just as a safety check (even though your query should enforce it).
    """
    dimers = {}
    skipped = 0

    for pdb, d in per_pdb.items():
        seqs = d["seqs"]
        if len(seqs) != 2:
            skipped += 1
            continue

        # Deterministic ordering by entity ID so results are stable
        ent_ids = sorted(seqs.keys(), key=lambda x: int(x) if x.isdigit() else x)
        s1, s2 = seqs[ent_ids[0]], seqs[ent_ids[1]]

        dimers[pdb] = (s1, s2, d["meta"], d["ligands"])

    return dimers, skipped


def cluster_exact_dimers(dimers):
    """
    Build exact sequence-pair clusters:
      clusters[key] = [pdb1, pdb2, ...]
    where key = order-independent (seqA, seqB).
    """
    clusters = defaultdict(list)
    for pdb, (s1, s2, meta, ligands) in dimers.items():
        clusters[dimer_key(s1, s2)].append(pdb)
    return clusters


def pick_representatives(clusters, dimers):
    """
    For each exact-sequence cluster, choose a representative using representative_score(meta).
    Returns:
      rep_for_cluster[key] = rep_pdb
    """
    rep_for_cluster = {}
    for key, members in clusters.items():
        rep = None
        best_score = (inf, inf)

        for pdb in members:
            meta = dimers[pdb][2]
            score = representative_score(meta)
            if score < best_score:
                best_score = score
                rep = pdb

        rep_for_cluster[key] = rep

    return rep_for_cluster


# =========================
# 5) REMOVE POINT-MUTANT DUPLICATES (AMONG REPS ONLY)
# =========================

def remove_point_mutant_duplicates(reps, dimers):
    """
    "Point mutant duplicate" definition (what I implement):
      - one partner sequence is identical
      - the other partner has same length and differs by exactly 1 residue (Hamming distance = 1)
      - no indels allowed (length mismatch => not considered)

    I decide which one to keep using the same representative_score(meta).
    """
    # Index representatives by shared partner sequence, so I only compare plausible pairs
    partner_index = defaultdict(list)  # partner_seq -> list of (rep_pdb, other_seq)
    for pdb in reps:
        a, b, meta, ligands = dimers[pdb]
        partner_index[b].append((pdb, a))
        partner_index[a].append((pdb, b))

    to_remove = set()
    decision_log = []

    for partner_seq, items in partner_index.items():
        for (pdb1, s1), (pdb2, s2) in combinations(items, 2):
            d = hamming(s1, s2)
            if d == 1:
                # Keep the one with the "better" score
                score1 = representative_score(dimers[pdb1][2])
                score2 = representative_score(dimers[pdb2][2])

                keep = pdb1 if score1 <= score2 else pdb2
                drop = pdb2 if keep == pdb1 else pdb1

                to_remove.add(drop)

                decision_log.append({
                    "keep_pdb": keep,
                    "drop_pdb": drop,
                    "keep_ligand_count": dimers[keep][2].get("ligand_count"),
                    "drop_ligand_count": dimers[drop][2].get("ligand_count"),
                    "keep_quality": quality_score(dimers[keep][2]),
                    "drop_quality": quality_score(dimers[drop][2]),
                    "partner_sequence_len": len(partner_seq),
                })

    final_reps = sorted(set(reps) - to_remove)
    return final_reps, decision_log


# =========================
# 6) WRITE OUTPUT FILES
# =========================

def write_primary_ids(pdb_ids, out_path):
    """Write a one-column CSV containing only the final PDB IDs."""
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["PDB ID"])
        for pdb in pdb_ids:
            w.writerow([pdb])

def write_exact_cluster_mapping(clusters, rep_for_cluster, dimers, out_path):
    """
    For each exact sequence-pair cluster, write:
      - cluster_id
      - representative_pdb
      - member_pdb
      - metadata (method/resolution)
      - ligand_count and ligand_names (joined)
      - cluster size

    This is the file that lets me "bring back" ligand-bound alternates later,
    without having to redo the whole search.
    """
    fieldnames = [
        "cluster_id",
        "representative_pdb",
        "member_pdb",
        "member_method",
        "member_em_res",
        "member_xray_res",
        "member_quality",
        "member_ligand_count",
        "member_ligand_names",
        "n_members_in_cluster",
    ]

    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()

        for key, members in clusters.items():
            cid = cluster_id_from_key(key)
            rep = rep_for_cluster[key]
            n = len(members)

            for pdb in sorted(members):
                s1, s2, meta, ligands = dimers[pdb]
                w.writerow({
                    "cluster_id": cid,
                    "representative_pdb": rep,
                    "member_pdb": pdb,
                    "member_method": meta.get("method", ""),
                    "member_em_res": meta.get("em_res", ""),
                    "member_xray_res": meta.get("xray_res", ""),
                    "member_quality": quality_score(meta) if quality_score(meta) != inf else "",
                    "member_ligand_count": meta.get("ligand_count", ""),
                    "member_ligand_names": "; ".join(sorted(ligands)) if ligands else "",
                    "n_members_in_cluster": n,
                })

def write_point_mutant_log(decisions, out_path):
    """Write a log of what I removed as point-mutant duplicates and why."""
    fieldnames = [
        "keep_pdb",
        "drop_pdb",
        "keep_ligand_count",
        "drop_ligand_count",
        "keep_quality",
        "drop_quality",
        "partner_sequence_len",
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in decisions:
            w.writerow(row)


# =========================
# 7) MAIN RUN
# =========================

def main():
    # (A) Combine chunks into one file
    input_files = combine_reports(INPUT_GLOB, COMBINED_CSV)
    print(f"Combined {len(input_files)} files into {COMBINED_CSV}")

    # (B) Parse the combined report
    per_pdb = parse_combined_report(COMBINED_CSV)

    # (C) Build (seqA, seqB) per PDB ID
    dimers, skipped = build_dimers(per_pdb)

    # (D) Exact clustering by sequence pair
    clusters = cluster_exact_dimers(dimers)

    # (E) Pick one representative per exact cluster
    rep_for_cluster = pick_representatives(clusters, dimers)
    exact_reps = sorted(set(rep_for_cluster.values()))

    # (F) Remove point-mutant duplicates among those reps
    final_primary, point_mutant_decisions = remove_point_mutant_duplicates(exact_reps, dimers)

    # (G) Write outputs
    write_primary_ids(final_primary, OUT_PRIMARY)
    write_exact_cluster_mapping(clusters, rep_for_cluster, dimers, OUT_CLUSTERS)
    write_point_mutant_log(point_mutant_decisions, OUT_POINT_MUT)

    # (H) Print a short summary so I can sanity-check counts
    print(f"Parsed PDB entries: {len(per_pdb)}")
    print(f"Entries with 2 entity sequences: {len(dimers)} (skipped {skipped})")
    print(f"Exact clusters (unique sequence-pairs): {len(clusters)}")
    print(f"Exact representatives: {len(exact_reps)}")
    print(f"Primary after point-mutant pruning: {len(final_primary)}")
    print(f"Wrote: {OUT_PRIMARY}, {OUT_CLUSTERS}, {OUT_POINT_MUT}")

if __name__ == "__main__":
    main()
