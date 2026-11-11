#!/usr/bin/env python3
"""Filter phaser_cis_var output and add Benjamini-Hochberg FDR."""

import argparse
import csv
from pathlib import Path


def bh(pvalues):
    n = len(pvalues)
    out = [1.0] * n
    previous = 1.0
    for rank, index in reversed(list(enumerate(sorted(range(n), key=pvalues.__getitem__), 1))):
        previous = min(previous, pvalues[index] * n / rank)
        out[index] = previous
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="phaser_cis_var.py TSV output")
    parser.add_argument("output", help="filtered output TSV")
    parser.add_argument("--min-hets", type=int, default=10)
    parser.add_argument("--p", type=float, default=0.05)
    parser.add_argument("--fdr", type=float, default=0.05)
    args = parser.parse_args()

    with open(args.input, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"var_het_n", "var_het_pval"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit(f"Input must contain: {', '.join(sorted(required))}")

    pvalues = [float(row["var_het_pval"]) for row in rows]
    fdr = bh(pvalues)
    for row, value in zip(rows, fdr):
        row["var_het_fdr"] = f"{value:.6g}"
        row["pass_ASE"] = str(
            int(row["var_het_n"]) >= args.min_hets
            and float(row["var_het_pval"]) < args.p
            and value < args.fdr
        ).upper()

    selected = [row for row in rows if row["pass_ASE"] == "TRUE"]
    fields = list(rows[0])
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    print(f"tested={len(rows)} retained={len(selected)} output={args.output}")


if __name__ == "__main__":
    main()
