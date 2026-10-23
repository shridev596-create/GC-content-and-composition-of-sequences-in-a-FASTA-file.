#!/usr/bin/env python3
"""fasta_stats.py - report basic properties of sequences in a FASTA file.

For every record it reports:
  * length
  * GC content (%), computed over unambiguous bases (A, C, G, T/U)
  * composition (count and % of each character)

Pure standard library. Works on any FASTA file (DNA or RNA).

Usage:
    python fasta_stats.py input.fasta
    python fasta_stats.py input.fasta --out report.txt
    python fasta_stats.py input.fasta --csv report.csv
"""
import argparse
import csv
import sys
from collections import Counter


def parse_fasta(handle):
    """Yield (header, sequence) tuples from an open FASTA file handle.

    - Header is the text after '>' on the record line.
    - Sequence lines are joined, whitespace stripped, upper-cased.
    - Blank lines and ';' comment lines are ignored.
    """
    header = None
    chunks = []
    for raw in handle:
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks).upper()
            header = line[1:].strip()
            chunks = []
        else:
            if header is None:
                raise ValueError("File does not look like FASTA: sequence data before first '>' header")
            chunks.append(line.replace(" ", ""))
    if header is not None:
        yield header, "".join(chunks).upper()


def gc_content(seq):
    """GC% over A/C/G/T/U bases only (N and other symbols are excluded)."""
    counts = Counter(seq)
    gc = counts["G"] + counts["C"]
    total = gc + counts["A"] + counts["T"] + counts["U"]
    return 100.0 * gc / total if total else 0.0


def analyse(seq):
    comp = Counter(seq)
    length = len(seq)
    return {
        "length": length,
        "gc": gc_content(seq),
        "composition": dict(sorted(comp.items())),
    }


def format_record(header, stats):
    lines = [f">{header}", f"  Length     : {stats['length']} bp", f"  GC content : {stats['gc']:.2f} %"]
    comp = stats["composition"]
    parts = []
    for base, n in comp.items():
        pct = 100.0 * n / stats["length"] if stats["length"] else 0.0
        parts.append(f"{base}={n} ({pct:.2f}%)")
    lines.append("  Composition: " + ", ".join(parts))
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description="Basic stats for sequences in a FASTA file")
    ap.add_argument("fasta", help="path to FASTA file")
    ap.add_argument("--out", help="also write the text report to this file")
    ap.add_argument("--csv", help="write per-record table to this CSV file")
    args = ap.parse_args()

    results = []
    with open(args.fasta) as fh:
        for header, seq in parse_fasta(fh):
            results.append((header, analyse(seq)))

    if not results:
        sys.exit("No FASTA records found.")

    report = [f"File: {args.fasta}", f"Records: {len(results)}", ""]
    for header, stats in results:
        report.append(format_record(header, stats))
        report.append("")

    total_len = sum(s["length"] for _, s in results)
    all_comp = Counter()
    for _, s in results:
        all_comp.update(s["composition"])
    gc = all_comp["G"] + all_comp["C"]
    acgt = gc + all_comp["A"] + all_comp["T"] + all_comp["U"]
    report.append("=== Overall ===")
    report.append(f"Total length : {total_len} bp")
    report.append(f"Overall GC   : {100.0 * gc / acgt if acgt else 0.0:.2f} %")
    text = "\n".join(report)
    print(text)

    if args.out:
        with open(args.out, "w") as f:
            f.write(text + "\n")

    if args.csv:
        chars = sorted(all_comp)
        with open(args.csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["id", "length", "gc_percent"] + chars)
            for header, s in results:
                w.writerow([header, s["length"], f"{s['gc']:.2f}"] + [s["composition"].get(c, 0) for c in chars])


if __name__ == "__main__":
    main()
