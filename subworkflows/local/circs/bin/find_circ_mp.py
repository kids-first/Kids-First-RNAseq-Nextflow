#!/usr/bin/env python3
import argparse
import codecs
import logging
import mmap
import os
import sys
import traceback
from collections import defaultdict
from itertools import zip_longest

import numpy
import pysam
from numpy import byte, frombuffer

import multiprocessing as mp
from types import SimpleNamespace
from typing import Any


_WORKER_OPTIONS = None
_WORKER_GENOME = None


class indexed_fasta:
    def __init__(self, fname: str, split_chrom: str = "", **_kwargs):
        self.fname = fname
        self.chrom_stats: dict[str, tuple[int, int, int, str, int]] = {}
        self.split_chrom = split_chrom

        ipath = fname + ".byo_index"
        if os.access(ipath, os.R_OK):
            self.load_index(ipath)
        else:
            self.index()
            self.store_index(ipath)

        self._f_handle = open(fname, "rb")  # keep alive for mmap lifetime
        self.mmap = mmap.mmap(self._f_handle.fileno(), 0, access=mmap.ACCESS_READ)

    def close(self) -> None:
        try:
            self.mmap.close()
        finally:
            self._f_handle.close()

    def index(self) -> None:
        logging.debug("# indexed_fasta.index(%r)", self.fname)

        chrom = "undef"
        chrom_ofs = 0
        size = 0
        ofs = 0

        stats_build: dict[str, list] = {}
        with open(self.fname, "rb") as f:
            for line_b in f:
                ofs += len(line_b)
                if line_b.startswith(b">"):
                    if size and chrom in stats_build:
                        stats_build[chrom].append(size)

                    chrom = line_b[1:].split()[0].decode("ascii", "replace").strip()
                    if self.split_chrom:
                        chrom = chrom.split(self.split_chrom)[0]

                    chrom_ofs = ofs
                    size = 0
                else:
                    if chrom not in stats_build:
                        lline = len(line_b)
                        stripped = line_b.rstrip(b"\r\n")
                        ldata = len(stripped)
                        nl_char = lline - ldata
                        skipchar = line_b[ldata:lline].decode("ascii", "replace")
                        stats_build[chrom] = [chrom_ofs, ldata, nl_char, skipchar]
                    size += len(line_b.rstrip(b"\r\n"))

        if size and chrom in stats_build:
            stats_build[chrom].append(size)

        self.chrom_stats = {}
        for chrom, vals in stats_build.items():
            chrom_ofs, ldata, skip, skipchar, size = vals
            self.chrom_stats[chrom] = (int(chrom_ofs), int(ldata), int(skip), str(skipchar), int(size))

    def store_index(self, ipath: str) -> None:
        logging.debug("# indexed_fasta.store_index(%r)", ipath)

        import stat
        import tempfile

        tmp = tempfile.NamedTemporaryFile(mode="w", dir=os.path.dirname(ipath) or ".", delete=False)
        try:
            for chrom in sorted(self.chrom_stats.keys()):
                ofs, ldata, skip, skipchar, size = self.chrom_stats[chrom]
                tmp.write(f"{chrom}\t{ofs}\t{ldata}\t{skip}\t{skipchar!r}\t{size}\n")
            tmp.flush()
            os.fsync(tmp.fileno())
        finally:
            tmp.close()

        os.chmod(tmp.name, stat.S_IROTH | stat.S_IRGRP | stat.S_IRUSR)
        os.replace(tmp.name, ipath)

    def load_index(self, ipath: str) -> None:
        logging.debug("# indexed_fasta.load_index(%r)", ipath)
        self.chrom_stats = {}
        with open(ipath, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                chrom, ofs, ldata, skip, skipchar, size = line.rstrip("\n").split("\t")
                decoded_skipchar = codecs.decode(skipchar[1:-1], "unicode_escape")
                self.chrom_stats[chrom] = (int(ofs), int(ldata), int(skip), decoded_skipchar, int(size))

    def get_data(self, chrom: str, start: int, end: int, sense: str) -> str:
        if not self.chrom_stats:
            self.index()

        ofs, ldata, skip, skip_char, size = self.chrom_stats[chrom]
        pad_start = 0
        pad_end = 0

        if start < 0:
            pad_start = -start
            start = 0
        if end > size:
            pad_end = end - size
            end = size

        l_start = start // ldata
        l_end = end // ldata
        ofs_start = l_start * skip + start + ofs
        ofs_end = l_end * skip + end + ofs

        chunk = self.mmap[ofs_start:ofs_end]
        s = chunk.decode("ascii", "replace").replace(skip_char, "")

        if pad_start or pad_end:
            s = ("N" * pad_start) + s + ("N" * pad_end)
        if sense == "-":
            s = rev_comp(s)
        return s


class Accessor:
    supports_write = False

    def __init__(self, _path: str, chrom: str, sense: str, sense_specific: bool = False, **_kwargs):
        if sense_specific:
            self.covered_strands = [chrom + sense]
        else:
            self.covered_strands = [chrom + "+", chrom + "-"]

    def get_data(self, chrom: str, start: int, end: int, sense: str, **_kwargs):
        return []

    def get_oriented(self, chrom: str, start: int, end: int, sense: str, **kwargs):
        data = self.get_data(chrom, start, end, sense, **kwargs)
        return data[::-1] if sense == "-" else data

    def flush(self) -> None:
        return


class Track:
    def __init__(
        self,
        path: str,
        accessor,
        sense_specific: bool = True,
        description: str = "unlabeled track",
        system: str = "hg18",
        dim: int = 1,
        auto_flush: bool = False,
        mode: str = "r",
        **kwargs,
    ):
        self.path = path
        self.mode = mode
        self.acc_cache: dict[str, Accessor] = {}
        self.accessor = accessor
        self.kwargs = kwargs
        self.sense_specific = bool(sense_specific)
        self.dim = int(dim)
        self.description = description
        self.auto_flush = auto_flush
        self.last_chrom = ""
        self.logger = logging.getLogger(f"Track({path!r})")
        self.system = system

        kwargs["sense_specific"] = self.sense_specific
        kwargs["mode"] = self.mode
        kwargs["system"] = self.system
        kwargs["description"] = self.description
        kwargs["dim"] = self.dim

    def load(self, chrom: str, sense: str) -> Accessor:
        if self.auto_flush and chrom != self.last_chrom:
            self.logger.debug("Seen new chromosome %s. Flushing accessor caches.", chrom)
            self.flush_all()
        self.last_chrom = chrom

        ID = chrom + sense
        if ID not in self.acc_cache:
            self.logger.debug("Cache miss for %s%s. creating new accessor", chrom, sense)
            acc = self.accessor(self.path, chrom, sense, **self.kwargs)

            if acc.covered_strands == "*":
                self.logger.debug("Accessor claims global coverage; registering for %s only (safe fallback).", ID)
                self.acc_cache[ID] = acc
            else:
                for sid in acc.covered_strands:
                    self.acc_cache[sid] = acc

        return self.acc_cache[ID]

    def flush_all(self) -> None:
        for a in set(self.acc_cache.values()):
            a.flush()
        self.acc_cache = {}

    def get(self, chrom: str, start: int, end: int, sense: str, **kwargs):
        return self.load(chrom, sense).get_data(chrom, start, end, sense, **kwargs)

    def get_oriented(self, chrom: str, start: int, end: int, sense: str, **kwargs):
        return self.load(chrom, sense).get_oriented(chrom, start, end, sense, **kwargs)


class GenomeAccessor(Accessor):
    def __init__(self, path: str, chrom: str, sense: str, system: str = "hg19", **kwargs):
        super().__init__(path, chrom, sense, system=system, **kwargs)
        logging.debug("# GenomeAccessor mmap: Loading genomic sequence for chromosome %s from %r", chrom, path)

        self.system = system
        try:
            self.data = indexed_fasta(os.path.join(path))
        except OSError:
            logging.warning("Could not access %r. Switching to dummy mode (only Ns)", path)
            self.data = None
            self.get_data = self.get_dummy  # type: ignore[method-assign]
            self.get_oriented = self.get_dummy  # type: ignore[method-assign]
            self.covered_strands = [chrom + "+", chrom + "-"]
        else:
            chroms = list(self.data.chrom_stats.keys())
            self.covered_strands = [c + "+" for c in chroms] + [c + "-" for c in chroms]

        self.get = self.get_oriented  # legacy alias; keep if other code expects it

    def get_data(self, chrom: str, start: int, end: int, sense: str):
        seq = self.data.get_data(chrom, start, end, "+")  # type: ignore[union-attr]
        if sense == "-":
            seq = complement(seq)  # matches the original script behavior
        return seq

    def get_dummy(self, _chrom: str, start: int, end: int, _sense: str):
        return "N" * int(end - start)


class Hit:
    def __init__(self, options, samples, minmapscore: int):
        self.options = options
        self.samples = samples
        self.minmapscore = minmapscore

        self.reads: list[str] = []
        self.readnames: list[str] = []
        self.uniq: set[tuple[str, str]] = set()
        self.mapquals_A: list[int] = []
        self.mapquals_B: list[int] = []
        self.uniq_bridges = 0
        self.tissues = defaultdict(int)
        self.edits: list[int] = []
        self.overlaps: list[int] = []
        self.n_hits: list[int] = []
        self.signal = "NNNN"
        self.strand_plus = 0
        self.strand_minus = 0
        self.strandmatch = "NA"

    def add(self, read: str, A, B, dist: int, ov: int, strandmatch: str, signal: str, n_hits: int) -> None:
        self.signal = signal
        self.strandmatch = strandmatch
        self.edits.append(dist)
        self.overlaps.append(ov)
        self.n_hits.append(n_hits)

        if A.pos > B.pos:
            A, B = B, A

        aopt = dict(A.tags)
        bopt = dict(B.tags)
        qA = aopt.get("AS") - aopt.get("XS", self.minmapscore)
        qB = bopt.get("AS") - bopt.get("XS", self.minmapscore)

        if qA and qB:
            self.uniq_bridges += 1

        self.mapquals_A.append(qA)
        self.mapquals_B.append(qB)

        if "__" in A.query_name:
            qname = A.query_name.split("__")[0][:-2]
        else:
            qname = B.query_name.split("__")[0][:-2]

        self.readnames.append(qname)

        if A.is_reverse:
            self.strand_minus += 1
            self.reads.append(rev_comp(read))
        else:
            self.strand_plus += 1
            self.reads.append(read)

        sample_name = self.options.name
        for (prefix, tiss) in self.samples:
            if qname.startswith(prefix):
                sample_name = tiss
                break

        self.tissues[sample_name] += 1
        self.uniq.add((read, sample_name))
        self.uniq.add((rev_comp(read), sample_name))

    def scores(self):
        n_reads = len(self.reads)
        n_uniq = len(self.uniq) // 2
        best_qual_A = sorted(self.mapquals_A, reverse=True)[0]
        best_qual_B = sorted(self.mapquals_B, reverse=True)[0]
        tissues = sorted(self.tissues.keys())
        tiss_counts = [str(self.tissues[k]) for k in tissues]
        return (
            n_reads,
            n_uniq,
            best_qual_A,
            best_qual_B,
            self.uniq_bridges,
            tissues,
            tiss_counts,
            min(self.edits),
            min(self.overlaps),
            min(self.n_hits),
            self.signal,
            self.strandmatch,
        )


COMPLEMENT = {
    "a": "t", "t": "a", "c": "g", "g": "c", "k": "m", "m": "k", "r": "y", "y": "r",
    "s": "s", "w": "w", "b": "v", "v": "b", "h": "d", "d": "h", "n": "n",
    "A": "T", "T": "A", "C": "G", "G": "C", "K": "M", "M": "K", "R": "Y", "Y": "R",
    "S": "S", "W": "W", "B": "V", "V": "B", "H": "D", "D": "H", "N": "N",
}


def complement(s: str) -> str:
    return "".join(COMPLEMENT[x] for x in s)


def rev_comp(seq: str) -> str:
    return complement(seq)[::-1]


def find_breakpoints(options, genome: Track, A, B, read: str, chrom: str):
    def strandmatch(ann: str, sense: str) -> str:
        if ann == sense:
            return "MATCH"
        if ann == "*" or len(ann) > 1:
            return "NA"
        return "MISMATCH"

    rnd = (lambda: float(numpy.random.random())) if options.randomize else (lambda: 0.0)

    L = len(read)
    hits = []
    eff_a = options.asize - options.margin
    internal = read[eff_a:-eff_a].upper()

    flank = L - 2 * eff_a + 2
    A_flank = genome.get(
        chrom,
        A.reference_end - options.margin,
        A.reference_end - options.margin + flank,
        "+",
    ).upper()
    B_flank = genome.get(
        chrom,
        B.reference_start - flank + options.margin,
        B.reference_start + options.margin,
        "+",
    ).upper()

    l = L - 2 * eff_a  # len(internal)
    if l < 0:
        return hits

    # Precompute mismatch prefix/suffix counts (no per-x spliced construction)
    I = frombuffer(internal.encode("ascii", "replace"), dtype=byte)
    Af = frombuffer(A_flank.encode("ascii", "replace"), dtype=byte)
    Bf = frombuffer(B_flank.encode("ascii", "replace"), dtype=byte)

    left_eq = (I[:l] == Af[:l])  # length l
    left_mis = numpy.empty(l + 1, dtype=numpy.int32)
    left_mis[0] = 0
    left_mis[1:] = numpy.arange(1, l + 1, dtype=numpy.int32) - numpy.cumsum(left_eq, dtype=numpy.int32)

    right_eq = (I[:l] == Bf[2 : 2 + l])  # I[i] vs Bf[i+2]
    right_mis = numpy.empty(l + 1, dtype=numpy.int32)
    right_mis[l] = 0
    right_mis[:l] = numpy.cumsum(~right_eq[::-1], dtype=numpy.int32)[::-1]

    # Candidate x positions:
    # - noncanonical: evaluate all x
    # - canonical: only x where A_flank[x:x+2] + B_flank[x:x+2] is GTAG or CTAC
    if options.noncanonical:
        x_candidates = range(l + 1)
    else:
        x_candidates = []
        max_x = l
        for x in range(max_x + 1):
            gt = A_flank[x : x + 2]
            ag = B_flank[x : x + 2]
            if gt == "GT" and ag == "AG":
                x_candidates.append(x)  # GTAG (+)
            elif gt == "CT" and ag == "AC":
                x_candidates.append(x)  # CTAC (-)

    for x in x_candidates:
        dist = int(left_mis[x] + right_mis[x])

        ov = 0
        if x < options.margin:
            ov = options.margin - x
        if l - x < options.margin:
            ov = max(ov, options.margin - (l - x))

        if dist > options.maxdist:
            continue

        gt = A_flank[x : x + 2]
        ag = B_flank[x : x + 2]
        gtag = gt + ag
        rc_gtag = rev_comp(gtag)

        start = B.reference_start + options.margin - l + x
        end = A.reference_end - options.margin + x + 1
        start, end = min(start, end), max(start, end)

        strand = "*"
        if options.stranded:
            strand = "-" if A.is_reverse else "+"

        if options.noncanonical:
            hits.append((dist, ov, strandmatch(strand, "+"), rnd(), chrom, start, end, gtag, "+"))
            hits.append((dist, ov, strandmatch(strand, "-"), rnd(), chrom, start, end, rc_gtag, "-"))
        else:
            if gtag == "GTAG":
                hits.append((dist, ov, strandmatch(strand, "+"), rnd(), chrom, start, end, "GTAG", "+"))
            elif gtag == "CTAC":
                hits.append((dist, ov, strandmatch(strand, "-"), rnd(), chrom, start, end, "GTAG", "-"))

    if len(hits) < 2:
        return hits

    hits = sorted(hits)
    best = hits[0]
    if options.strandpref:
        return [h for h in hits if (h[0] == best[0]) and (h[1] == best[1]) and (h[2] == best[2])]
    return [h for h in hits if (h[0] == best[0]) and (h[1] == best[1])]


def grouper(n: int, iterable, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def output(options, bedfile, readfile, N, cand: dict, prefix: str) -> None:
    n = 1
    for (chrom, start, end, sense), hit in cand.items():
        (
            n_reads, n_uniq, best_qual_A, best_qual_B, uniq_bridges,
            tissues, tiss_counts, min_edit, min_anchor_ov, n_hits, signal, strandmatch,
        ) = hit.scores()

        if options.halfunique:
            if (best_qual_A < options.min_uniq_qual) and (best_qual_B < options.min_uniq_qual):
                N["anchor_not_uniq"] += 1
                continue
        else:
            if (best_qual_A < options.min_uniq_qual) or (best_qual_B < options.min_uniq_qual):
                N["anchor_not_uniq"] += 1
                continue

        if (uniq_bridges == 0) and (not options.report_nobridges):
            N["no_uniq_bridges"] += 1
            continue

        name = f"{options.prefix}{prefix}_{n:06d}"
        n += 1

        for r_seq, ori_name in zip(hit.reads, hit.readnames):
            print(f">{name} {ori_name}", file=readfile)
            print(r_seq, file=readfile)

        categories = []
        if signal == "GTAG":
            categories.append("CANONICAL")
        if strandmatch == "MATCH":
            categories.append("STRANDMATCH")
        if best_qual_A > 0 and best_qual_B > 0 and uniq_bridges > 0:
            categories.append("ANCHOR_UNIQUE")
        if uniq_bridges == 0:
            categories.append("NO_UNIQ_BRIDGES")
        if n_hits == 1:
            categories.append("UNAMBIGUOUS_BP")
        if min_anchor_ov == 0 and min_edit == 0:
            categories.append("PERFECT_EXT")
        elif min_anchor_ov <= 1 and min_edit <= 1:
            categories.append("GOOD_EXT")
        elif min_anchor_ov <= 2 and min_edit <= 2:
            categories.append("OK_EXT")
        if not categories:
            categories.append("DUBIOUS")

        categories.append("CIRCULAR" if prefix == "circ" else "LINEAR")

        bed = [
            chrom,
            start - 1,
            end,
            name,
            n_reads,
            sense,
            n_uniq,
            uniq_bridges,
            best_qual_A,
            best_qual_B,
            ",".join(tissues),
            ",".join(tiss_counts),
            min_edit,
            min_anchor_ov,
            n_hits,
            signal,
            strandmatch,
            ",".join(sorted(categories)),
        ]
        print("\t".join(map(str, bed)), file=bedfile)


def build_argparser() -> argparse.ArgumentParser:
    usage = "bowtie2 [mapping options] anchors.fastq.gz | %(prog)s [options] > candidates.bed"
    p = argparse.ArgumentParser(usage=usage)

    p.add_argument("-v", "--version", action="store_true", help="get version information")
    p.add_argument("-S", "--system", default="", help="model system database (optional! Requires byo library.)")
    p.add_argument("-G", "--genome", default="", help="path to genome (folder with chr*.fa or one multi-chrom FASTA)")
    p.add_argument("-n", "--name", default="unknown", help="tissue/sample name (default: unknown)")
    p.add_argument("-p", "--prefix", default="", help="prefix to prepend to each junction name")
    p.add_argument("-q", "--min_uniq_qual", type=int, default=2, help="minimal uniqueness for anchor alignments")
    p.add_argument("-a", "--anchor", dest="asize", type=int, default=20, help="anchor size")
    p.add_argument("-m", "--margin", type=int, default=2, help="max nts breakpoint allowed inside an anchor")
    p.add_argument("-d", "--maxdist", type=int, default=2, help="max mismatches allowed in anchor extensions")

    p.add_argument("--noncanonical", action="store_true", help="relax the GU/AG constraint")
    p.add_argument("--randomize", action="store_true", help="select randomly from tied hits")
    p.add_argument("--allhits", action="store_true", help="in case of ambiguities, report each hit")
    p.add_argument("--stranded", action="store_true", help="use if the reads are stranded")
    p.add_argument("--strandpref", action="store_true", help="prefer splice sites matching annotated transcription")
    p.add_argument("--halfunique", action="store_true", help="report junctions where only one anchor aligns uniquely")
    p.add_argument("--report_nobridges", action="store_true", help="also report junctions lacking uniq bridges")

    p.add_argument("--output_prefix", required=True, help="prefix for all outputs (required)")
    p.add_argument("--write-bam", action="store_true", help="write candidate anchor alignments BAM")
    p.add_argument("--write-stats", action="store_true", help="write run statistics file")
    p.add_argument("-r", "--reads2samples", default="", help="TSV: read-name prefix -> sample ID mapping")

    p.add_argument("--procs", type=int, default=0, help="number of worker processes for breakpoint search (0=off)")
    p.add_argument("--batch", type=int, default=512, help="batch size (number of breakpoint tasks per dispatch)")
    return p


def _init_worker(options_dict: dict[str, Any]) -> None:
    """
    Per-process initializer. Builds genome Track once in each worker.
    """
    global _WORKER_OPTIONS, _WORKER_GENOME

    _WORKER_OPTIONS = SimpleNamespace(**options_dict)

    genome = None
    if _WORKER_OPTIONS.system:
        import importlib
        system = importlib.import_module(f"byo.systems.{_WORKER_OPTIONS.system}")
        genome = system.genome

    if _WORKER_OPTIONS.genome:
        genome = Track(_WORKER_OPTIONS.genome, accessor=GenomeAccessor)

    if genome is None:
        raise RuntimeError("Worker could not initialize genome; specify -S or -G.")

    _WORKER_GENOME = genome


def _anchor_to_info(A) -> dict[str, Any]:
    """
    Extract only what find_breakpoints() needs from pysam.AlignedSegment.
    Must remain picklable.
    """
    return {
        "reference_start": int(A.reference_start),
        "reference_end": int(A.reference_end),
        "is_reverse": bool(A.is_reverse),
        "tags": list(A.tags),
        "query_name": str(A.query_name),
    }


def _find_breakpoints_task(payload: tuple[dict[str, Any], dict[str, Any], str, str]):
    """
    Worker-side wrapper. Reconstructs minimal A/B objects and calls find_breakpoints().
    """
    global _WORKER_OPTIONS, _WORKER_GENOME
    A_info, B_info, read, chrom = payload

    A = SimpleNamespace(**A_info)
    B = SimpleNamespace(**B_info)

    return find_breakpoints(_WORKER_OPTIONS, _WORKER_GENOME, A, B, read, chrom)


def main() -> int:
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    parser = build_argparser()
    options = parser.parse_args()

    if options.version:
        print("find_circ.py version 1.2\n\n(c) Marvin Jens 2012-2015.\nCheck http://www.circbase.org for more information.")
        return 0

    genome = None
    if options.system:
        import importlib
        system = importlib.import_module(f"byo.systems.{options.system}")
        genome = system.genome

    if options.genome:
        genome = Track(options.genome, accessor=GenomeAccessor)

    if genome is None:
        logging.error("need to specify either model system database (-S) or genome (-G).")
        return 1

    if options.reads2samples:
        with open(options.reads2samples, "r", encoding="utf-8", errors="replace") as f:
            samples = [line.rstrip("\n").split("\t") for line in f]
    else:
        samples = []
    samples.append(("", options.name))

    minmapscore = options.asize * (-2)
    circs = defaultdict(lambda: Hit(options, samples, minmapscore))
    splices = defaultdict(lambda: Hit(options, samples, minmapscore))
    N = defaultdict(int)

    sam = pysam.AlignmentFile("-", "r")
    bam_out = pysam.AlignmentFile(f"{options.output_prefix}.anchors.bam", "wb", template=sam) if options.write_bam else None

    ctx = mp.get_context("spawn")
    pool = None
    if options.procs and options.procs > 0:
        # Pass only simple types into initializer
        options_dict = vars(options).copy()
        pool = ctx.Pool(processes=options.procs, initializer=_init_worker, initargs=(options_dict,))

    pending_payloads: list[tuple[dict[str, Any], dict[str, Any], str, str]] = []
    pending_meta: list[tuple[str, str, Any, Any, str]] = []
    # meta tuple: (kind, chrom, A, B, read) where kind in {"circ","splice"}

    def flush_pending() -> None:
        nonlocal pending_payloads, pending_meta

        if not pending_payloads:
            return

        if pool is None:
            results = list(map(_find_breakpoints_task, pending_payloads))
        else:
            results = pool.map(_find_breakpoints_task, pending_payloads)

        for bp, meta in zip(results, pending_meta):
            kind, chrom, A, B, read = meta

            if kind == "circ":
                if not bp:
                    N["circ_no_bp"] += 1
                else:
                    N["circ_reads"] += 1

                n_hits = len(bp)
                if bp and not options.allhits:
                    bp = [bp[0]]

                for h in bp:
                    dist, ov, strandmatch, _rnd, _chrom, start, end, signal, sense = h
                    key = (chrom, start + 1, end - 1, sense)
                    circs[key].add(read, A, B, dist, ov, strandmatch, signal, n_hits)

            else:  # "splice"
                if not bp:
                    N["splice_no_bp"] += 1
                else:
                    N["spliced_reads"] += 1

                n_hits = len(bp)
                if bp and not options.allhits:
                    bp = [bp[0]]

                for h in bp:
                    dist, ov, strandmatch, _rnd, _chrom, start, end, signal, sense = h
                    key = (chrom, start, end, sense)
                    splices[key].add(read, A, B, dist, ov, strandmatch, signal, n_hits)

        pending_payloads = []
        pending_meta = []

    pair_num = 0
    try:
        for pair_num, (A, B) in enumerate(grouper(2, sam)):
            if A is None or B is None:
                break

            N["total"] += 1
            if A.is_unmapped or B.is_unmapped:
                N["unmapped"] += 1
                continue
            if A.reference_id != B.reference_id:
                N["other_chrom"] += 1
                continue
            if A.is_reverse != B.is_reverse:
                N["other_strand"] += 1
                continue

            dist = B.reference_start - A.reference_start
            if numpy.abs(dist) < options.asize:
                N["overlapping_anchors"] += 1
                continue

            if bam_out:
                bam_out.write(A)
                bam_out.write(B)

            chrom = sam.get_reference_name(A.reference_id)

            if (A.is_reverse and dist > 0) or ((not A.is_reverse) and dist < 0):
                read = A.query_name.split("__")[1]
                if A.is_reverse:
                    A, B = B, A
                    read = rev_comp(read)

                A_info = _anchor_to_info(A)
                B_info = _anchor_to_info(B)

                pending_payloads.append((A_info, B_info, read, chrom))
                pending_meta.append(("circ", chrom, A, B, read))

                if len(pending_payloads) >= options.batch:
                    flush_pending()

            elif (A.is_reverse and dist < 0) or ((not A.is_reverse) and dist > 0):
                read = A.query_name.split("__")[1]
                if A.is_reverse:
                    A, B = B, A
                    read = rev_comp(read)

                A_info = _anchor_to_info(A)
                B_info = _anchor_to_info(B)

                pending_payloads.append((A_info, B_info, read, chrom))
                pending_meta.append(("splice", chrom, A, B, read))

                if len(pending_payloads) >= options.batch:
                    flush_pending()

            else:
                N["fallout"] += 1
                logging.warning("unhandled read: A=%r B=%r", A, B)

    except KeyboardInterrupt:
        sam_line = pair_num * 2
        fastq_line = pair_num * 8
        logging.warning(
            "KeyboardInterrupt by user while processing alignment pair %d, on input starting at SAM line %d, FASTQ line %d",
            pair_num,
            sam_line,
            fastq_line,
        )
    except Exception:
        sam_line = pair_num * 2
        fastq_line = pair_num * 8
        logging.error(
            "Unhandled exception raised while processing alignment pair %d, on input starting at SAM line %d, FASTQ line %d",
            pair_num,
            sam_line,
            fastq_line,
        )
        traceback.print_exc(file=sys.stderr)
        return 1
    finally:
        if bam_out:
            bam_out.close()
        sam.close()
        try:
            flush_pending()
        except Exception:
            # If flushing fails during teardown, still proceed with closing resources
            logging.error("Failed to flush pending breakpoint tasks during shutdown.")
            traceback.print_exc(file=sys.stderr)
        if pool is not None:
            pool.close()
            pool.join()

    if options.write_stats:
        with open(f"{options.output_prefix}.stats.txt", "w", encoding="utf-8") as stats:
            print(str(dict(N)), file=stats)

    with open(f"{options.output_prefix}.sites.reads", "w", encoding="utf-8") as readfile, open(f"{options.output_prefix}.sites.bed", "w", encoding="utf-8") as bedfile:
        header = [
            "chrom", "start", "end", "name", "n_reads", "strand", "n_uniq", "uniq_bridges",
            "best_qual_left", "best_qual_right", "tissues", "tiss_counts", "edits",
            "anchor_overlap", "breakpoints", "signal", "strandmatch", "category",
        ]
        print("#", "\t".join(header), file=bedfile)
        output(options, bedfile, readfile, N, circs, "circ")
        output(options, bedfile, readfile, N, splices, "norm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
