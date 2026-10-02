#!/usr/bin/env python3
"""Identify which 10x SI-GA sample-index wells are present in a sequencing run.

Reads the index read (I1) FASTQ from a small, undemultiplexed test conversion,
counts 8-mers, and assigns each to a SI-GA well (<=1 mismatch to one of the
well's four oligos). Compares the result with the indices expected in
samples.csv and writes a per-run table.

usage: detect_10x_index.py --fastq I1.fastq.gz --table SI-GA_single_index.csv \
                           --samples samples.csv --run P165 --out detected_P165.csv
"""
import argparse
import csv
import gzip
import sys
from collections import Counter


def hamming1_neighbors(seq):
    out = {seq}
    for i, b in enumerate(seq):
        for alt in "ACGTN":
            if alt != b:
                out.add(seq[:i] + alt + seq[i + 1:])
    return out


def load_table(path):
    lookup, wells = {}, []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            well = row["index_name"]
            wells.append(well)
            for k in ("oligo1", "oligo2", "oligo3", "oligo4"):
                for s in hamming1_neighbors(row[k]):
                    # exact matches win over 1-mismatch matches
                    if s not in lookup or s == row[k]:
                        lookup[s] = well
    return lookup, wells


def count_index_reads(path, max_reads):
    counts = Counter()
    opener = gzip.open if path.endswith(".gz") else open
    n = 0
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                counts[line.strip()[:8]] += 1
                n += 1
                if n >= max_reads:
                    break
    return counts, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-reads", type=int, default=2_000_000)
    ap.add_argument("--min-fraction", type=float, default=0.005)
    a = ap.parse_args()

    lookup, _ = load_table(a.table)
    counts, n = count_index_reads(a.fastq, a.max_reads)
    if n == 0:
        sys.exit(f"[detect] no reads in {a.fastq}")

    per_well = Counter()
    for seq, c in counts.items():
        w = lookup.get(seq)
        if w:
            per_well[w] += c
    assigned = sum(per_well.values())

    expected = {}
    with open(a.samples) as fh:
        for row in csv.DictReader(fh):
            if row["run"].strip() == a.run:
                expected[row["index"].strip()] = row["sample_id"].strip()

    print(f"[detect] {a.run}: {n:,} index reads, {assigned / n:.1%} match a SI-GA well")
    rows = []
    for well, c in per_well.most_common():
        frac = c / n
        if frac < a.min_fraction and well not in expected:
            continue
        status = "expected" if well in expected else "UNEXPECTED"
        rows.append((well, expected.get(well, ""), c, frac, status))
    for well, sid in expected.items():
        if well not in per_well:
            rows.append((well, sid, 0, 0.0, "expected"))

    ok = True
    print(f"  {'well':<10}{'sample':<10}{'reads':>10}{'fraction':>10}  status")
    for well, sid, c, frac, status in rows:
        if status == "expected" and frac < a.min_fraction:
            status = "MISSING (expected but absent)"
            ok = False
        elif status == "UNEXPECTED":
            ok = False
        print(f"  {well:<10}{sid:<10}{c:>10,}{frac:>10.2%}  {status}")

    unk = Counter({s: c for s, c in counts.items() if s not in lookup})
    if unk and sum(unk.values()) / n > 0.2:
        print("  top unassigned I1 sequences (dual-index or non-10x library?):")
        for s, c in unk.most_common(8):
            print(f"    {s}  {c / n:.2%}")

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run", "index", "sample_id", "reads", "fraction", "status"])
        for well, sid, c, frac, status in rows:
            present = frac >= a.min_fraction
            w.writerow([a.run, well, sid, c, f"{frac:.5f}",
                        "present" if present else "absent"])
    if not ok:
        print(f"[detect] WARNING: {a.run} does not match samples.csv - check before step 03")


if __name__ == "__main__":
    main()
