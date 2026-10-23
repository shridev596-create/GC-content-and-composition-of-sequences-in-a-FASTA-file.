# FASTA Sequence Stats (Task 01 — Sequence handling)

A small, dependency-free Python script that reads a FASTA file and reports basic properties of every sequence in it.

## What it reports

For each record:
- **Length** (bp)
- **GC content** (%)
- **Composition** (count and % of every character, e.g. A, C, G, T, N)

Plus an overall summary (total length and overall GC) across all records.

## How to run

Requires Python 3.7+ only (standard library, nothing to install).

```bash
python fasta_stats.py sample.fasta
python fasta_stats.py your_file.fasta --out report.txt --csv report.csv
```

- `--out` saves the text report to a file
- `--csv` saves a per-record table (id, length, GC%, count of each character)

It works on any FASTA file: multi-line sequences, multiple records, blank lines, lowercase bases.

## Example output (`sample.fasta`, a tiny toy file)

See `sample_output.txt` and `sample_output.csv`.

```
>seq1 toy example, 40 bp
  Length     : 40 bp
  GC content : 50.00 %
  Composition: A=11 (27.50%), C=9 (22.50%), G=11 (27.50%), T=9 (22.50%)
```

## Decisions I made

- **GC content is calculated over A/C/G/T/U only.** Ambiguous bases such as `N` are counted in composition but left out of the GC denominator, so a few unknown bases don't drag the GC% down. Length still includes them.
- Sequences are upper-cased so `acgt` and `ACGT` are treated the same.
- FASTA is parsed with a script (no manual counting): headers start with `>`, sequence lines are joined until the next header.
- A file with no `>` header or no records stops with a clear error instead of printing nonsense.

## Running on a public dataset

Download any public FASTA (for example a genome or gene record from NCBI: open the record, then *Send to → File → FASTA*), then:

```bash
python fasta_stats.py downloaded.fasta --out output_public.txt
```

## Files

- `fasta_stats.py` — the script
- `sample.fasta` — tiny toy input to test quickly
- `sample_output.txt`, `sample_output.csv` — output for the toy input
