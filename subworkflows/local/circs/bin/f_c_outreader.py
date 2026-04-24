#!/usr/bin/env python3
# parsing the find_circ output : add coordinates in chr:start-end format, put out only important .bed columns

import argparse
import re
import sys


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-i", "--infile", required=True, help="Input file")
    p.add_argument("-o", "--outfile", default="out.bed", help="Output file")
    p.add_argument("--ignore_alt", "--alt", type=int, default=0, dest="ignore_alt", help="Include entries from alt contigs (e.g., chrUn_gl000220) if set to 0; exclude them if set to 1")
    p.add_argument("-s", "--strict", type=int, default=0, help="Only include entries with bestQA and bestQB >= 40")

    # Column indexes (0-based)
    p.add_argument("--idx-chr", type=int, default=0, help="Column index for chromosome")
    p.add_argument("--idx-beg", type=int, default=1, help="Column index for start position")
    p.add_argument("--idx-end", type=int, default=2, help="Column index for end position")

    p.add_argument("--idx-long-id", type=int, default=3, help="Column index for long ID")
    p.add_argument("--idx-strand", type=int, default=5, help="Column index for strand")
    p.add_argument("--idx-uniques", type=int, default=6, help="Column index for unique counts")
    p.add_argument("--idx-bestqa", type=int, default=7, help="Column index for bestQA")
    p.add_argument("--idx-bestqb", type=int, default=8, help="Column index for bestQB")

    # Scan for RefSeq starting at this column (0-based), inclusive
    p.add_argument("--idx-refseq-start", type=int, default=19, help="Column index for RefSeq start")

    return p.parse_args()

def pick_refseqid_from_parts(parts: list[str], start_index: int = 19) -> str:
    # Scan columns 19+ (0-based) and pick the first non-empty value that starts with "N"
    # (e.g., NM_, NR_, NP_, NG_).
    for v in parts[start_index:]:
        v = v.strip()
        if v.startswith("N"):
            return v
    raise ValueError(f"No valid RefseqID found in columns {start_index}+ for line: {'\t'.join(parts)}")

def main() -> int:
    args = parse_args()
    include_alts = bool(args.ignore_alt)

    with open(args.infile, "r", encoding="utf-8", errors="replace") as fin, open(
        args.outfile, "w", encoding="utf-8"
    ) as fout:
        print("\t".join(["coordinates", "strand", "sampleid", "unique_counts", "score", "score", "RefseqID"]), file=fout)

        for raw_line in fin:
            line = raw_line.rstrip("\n")
            if not line:
                continue

            parts = re.split(r"\t+", line)

            need_max = max(
                args.idx_chr,
                args.idx_beg,
                args.idx_end,
                args.idx_long_id,
                args.idx_strand,
                args.idx_uniques,
                args.idx_bestqa,
                args.idx_bestqb,
                args.idx_refseq_start,
            )
            if len(parts) <= need_max:
                print(f"Skipping line with insufficient columns (need at least {need_max + 1}): {line}", file=sys.stderr)
                continue

            chr = parts[args.idx_chr]
            beg = parts[args.idx_beg]
            end = parts[args.idx_end]
            ccord = f"{chr}:{beg}-{end}"

            # Strip "run_" and "circ_123456789" pattern from long_id to get a cleaner sample ID
            long_id = re.sub(r"run_", "", parts[args.idx_long_id], flags=re.IGNORECASE)
            long_id = re.sub(r"circ\_*.[0-9]{1,20}", "", long_id, flags=re.IGNORECASE)

            strand = parts[args.idx_strand]
            uniques = parts[args.idx_uniques]

            bestqa = parts[args.idx_bestqa]
            bestqb = parts[args.idx_bestqb]

            refseqid = pick_refseqid_from_parts(parts, start_index=args.idx_refseq_start)

            if args.strict and (bestqa < 40 or bestqb < 40):
                print(f"Skipping entry with bestQA={bestqa} and bestQB={bestqb} (below strict threshold of 40)", file=sys.stderr)
                continue

            if include_alts or not ccord.startswith("chrUn_gl"):
                print("\t".join([ccord, strand, long_id, uniques, bestqa, bestqb, refseqid]), file=fout)

if __name__ == "__main__":
    main()
