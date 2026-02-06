#!/usr/bin/env python
from collections import namedtuple, defaultdict
from itertools import izip_longest
from logging import debug, warning, getLogger
import multiprocessing as mp
import mmap
import numpy
from optparse import OptionParser
import os
import pysam
import sys
import traceback
import time
import json

# Worker-local globals (populated in worker_init)
WORKER_GENOME = None

COMPLEMENT = {
    'a': 't',
    't': 'a',
    'c': 'g',
    'g': 'c',
    'k': 'm',
    'm': 'k',
    'r': 'y',
    'y': 'r',
    's': 's',
    'w': 'w',
    'b': 'v',
    'v': 'b',
    'h': 'd',
    'd': 'h',
    'n': 'n',
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'K': 'M',
    'M': 'K',
    'R': 'Y',
    'Y': 'R',
    'S': 'S',
    'W': 'W',
    'B': 'V',
    'V': 'B',
    'H': 'D',
    'D': 'H',
    'N': 'N',
}

def complement(s):
    return "".join([COMPLEMENT[x] for x in s])

def rev_comp(seq):
    return complement(seq)[::-1]

def to_bool(obj):
    if obj == "False":
        return False
    else:
        return bool(obj)

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def mismatches_bytes_early(a_mv, b_mv, maxdist_local):
    # a_mv: memoryview/bytearray, b_mv: memoryview/bytes-like
    diff = 0
    ln = len(b_mv)
    # local variables for speed
    a_local = a_mv
    b_local = b_mv
    for i in xrange(ln):
        if a_local[i] != b_local[i]:
            diff += 1
            if diff > maxdist_local:
                return diff
    return diff

def find_breakpoints(A, B, read, genome, chrom, asize, margin, maxdist):
    L = len(read)
    hits = []
    eff_a = asize - margin
    internal = read[eff_a:-eff_a].upper()
    internal_b = internal.encode('ascii')
    flank = L - 2*eff_a + 2

    A_flank = genome.get(chrom, A.aend - margin, A.aend - margin + flank, '+').upper()
    B_flank = genome.get(chrom, B.pos - flank + margin, B.pos + margin, '+').upper()
    A_b = A_flank.encode('ascii')
    B_b = B_flank.encode('ascii')

    l = len(internal_b)
    spliced = bytearray(l)
    spliced_mv = memoryview(spliced)
    internal_mv = memoryview(internal_b)

    # local copies to avoid global lookups in hot loop
    mismatches_fn = mismatches_bytes_early
    maxdist_local = maxdist
    A_b_local = A_b
    B_b_local = B_b
    spliced_mv_local = spliced_mv
    internal_mv_local = internal_mv
    xrange_local = xrange

    for x in xrange_local(l + 1):
        # copy prefix from A_b and suffix from B_b in-place
        if x:
            spliced_mv_local[:x] = A_b_local[:x]
        suff_len = l - x
        if suff_len > 0:
            spliced_mv_local[x:] = B_b_local[x + 2: x + 2 + suff_len]

        # full (early-exit) mismatch count
        dist = mismatches_fn(spliced_mv_local, internal_mv_local, maxdist_local)

        ov = 0
        if x < margin:
            ov = margin - x
        if l - x < margin:
            ov = margin - (l - x)

        if dist <= maxdist_local:
            gt = A_flank[x:x+2]
            ag = B_flank[x:x+2]
            start, end = B.pos + margin - l + x, A.aend - margin + x + 1
            start, end = min(start, end), max(start, end)

            if gt == 'GT' and ag == 'AG':
                hits.append((dist, ov, chrom, start, end, '+'))
            elif gt == 'CT' and ag == 'AC':
                hits.append((dist, ov, chrom, start, end, '-'))

    if len(hits) < 2:
        return hits

    # Keep only ties on (edit distance, anchor overlap)
    hits = sorted(hits)
    best = hits[0]
    return [h for h in hits if (h[0] == best[0]) and (h[1] == best[1])]

def output(cand, prefix, suffix, wiggle, min_uniq_qual):
    print "#", "\t".join(['chrom','start','end','name','n_reads','strand','n_uniq','best_qual_A','best_qual_B','spliced_at_begin','spliced_at_end','tissues','tiss_counts','edits','anchor_overlap','breakpoints'])
    n = 1
    for c, hit in cand.items():
        chrom, start, end, sense = c
        n_reads, n_uniq, best_qual_A, best_qual_B, spliced_at_begin, spliced_at_end, tissues, tiss_counts, min_edit, min_anchor_ov, n_hits = hit.scores(chrom, start, end, sense, wiggle)

        if best_qual_B < min_uniq_qual:
            N['anchor_not_uniq'] += 1
            continue

        name = "%s%s_%06d" % (prefix, suffix, n)
        n += 1
        for r in hit.reads:
            sys.stderr.write("%s\t%s\n" % (name, r))

        bed = [chrom, start - 1, end, name, n_reads, sense, n_uniq, best_qual_A, best_qual_B, spliced_at_begin, spliced_at_end, ",".join(tissues), ",".join(tiss_counts), min_edit, min_anchor_ov, n_hits]
        print "\t".join([str(b) for b in bed])

class mmap_fasta(object):
    def __init__(self, fname):
        f = file(fname)
        header = f.readline()
        row = f.readline()

        self.ofs = len(header)
        self.lline = len(row)
        self.ldata = len(row.strip())
        self.skip = self.lline - self.ldata
        self.skip_char = row[self.ldata:]
        self.mmap = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)

    def __getslice__(self, start, end):
        l_start = start / self.ldata
        l_end = end / self.ldata
        ofs_start = l_start * self.skip + start + self.ofs
        ofs_end = l_end * self.skip + end + self.ofs

        s = self.mmap[ofs_start:ofs_end].replace(self.skip_char, "")
        L = end - start
        if len(s) == L:
            return s
        else:
            return s + "N" * (L - len(s))
        return

class Accessor(object):
    supports_write = False

    def get_data(self, start, end, sense):
        return []

    def get_oriented(self, start, end, sense):
        data = self.get_data(start, end, sense)
        if sense == "-":
            return data[::-1]
        else:
            return data

    def get_sum(self, start, end, sense):
        return self.get_data(start, end, sense).sum()

    def flush(self):
        pass

class Track(object):
    """
    Abstraction of chromosome-wide addressable data like sequences, coverage, scores etc.
    """

    def __init__(self, path, accessor, sense_specific=True, description="unlabeled track", system="hg19", dim=1, auto_flush=False, mode="r", **kwargs):
        self.path = path
        self.mode = mode
        self.acc_cache = {}
        self.accessor = accessor
        self.kwargs = kwargs
        self.sense_specific = to_bool(sense_specific)
        self.dim = int(dim)
        self.description = description
        self.auto_flush = auto_flush
        self.last_chrom = ""
        self.logger = getLogger("Track('%s')" % path)
        self.system = system

        self.logger.debug("Track(auto_flush=%s)" % (str(auto_flush)))
        kwargs['sense_specific'] = self.sense_specific
        kwargs['mode'] = self.mode
        kwargs['system'] = self.system
        kwargs['description'] = self.description
        kwargs['system'] = self.system
        kwargs['dim'] = self.dim

    def load(self, chrom, sense):
        if self.auto_flush and chrom != self.last_chrom:
            self.logger.debug("Seen new chromosome %s. Flushing accessor caches." % chrom)
            self.flush_all()
        self.last_chrom = chrom

        ID = self.get_identifier(chrom, sense)
        if not ID in self.acc_cache:
            self.logger.debug("Cache miss for %s%s. creating new accessor" % (chrom, sense))
            self.acc_cache[ID] = self.accessor(self.path, chrom, sense, **(self.kwargs))

        return self.acc_cache[ID]

    def save(self):
        if not "w" in self.mode:
            self.logger.warning("save() called on a read-only opened track. Ignored!")
            return

        if not self.accessor.supports_write:
            self.logger.warning("save() called on a track with only read-access supporting accessors. Ignored!")
            return

        self.logger.debug("save(): writing '%s'" % self.path)

        def to_str(obj):
            return getattr(obj, "__name__", str(obj))

        kwarg_str = "\n".join(["%s=%s" % (k, to_str(self.kwargs[k])) for k in sorted(self.kwargs.keys()) if k != "mode"])
        file(os.path.join(self.path, "track.rc"), "w+").write(trackrc % dict(accessor=self.accessor.__name__, kwargs=kwarg_str))
        self.flush_all()

    def __del__(self):
        if "w" in self.mode:
            self.save()

    def flush(self, chrom, sense):
        ID = self.get_identifier(chrom, sense)
        if ID in self.acc_cache:
            self.logger.warning("Flushing %s%s" % (chrom, sense))
            del self.acc_cache[ID]

    def flush_all(self):
        for a in self.acc_cache.values():
            a.flush()
        self.acc_cache = {}

    def get(self, chrom, start, end, sense):
        acc = self.load(chrom, sense)
        return acc.get_data(start, end, sense)

    def get_oriented(self, chrom, start, end, sense):
        acc = self.load(chrom, sense)
        return acc.get_oriented(start, end, sense)

    def get_sum(self, chrom, start, end, sense):
        acc = self.load(chrom, sense)
        return acc.get_sum(start, end, sense)

    def get_identifier(self, chrom, sense):
        if self.sense_specific:
            return chrom + sense
        else:
            return chrom

class GenomeAccessor(Accessor):
    def __init__(self, path, chrom, sense, **kwargs):
        debug("# GenomeAccessor mmap: Loading genomic sequence for chromosome %s from '%s'" % (chrom, path))

        fname = os.path.join(path, chrom + ".fa")
        try:
            self.data = mmap_fasta(fname)
        except IOError:
            warning("Could not access '%s'. Switching to dummy mode (only Ns)" % fname)
            self.get_data = self.get_dummy
            self.get_oriented = self.get_dummy

        self.get = self.get_oriented

    def get_data(self, start, end, sense):
        if start < 0 or end < 0:
            return self.get_dummy(start, end, sense)
        # UCSC convention: start with 1, end is inclusive
        if sense == "+":
            return self.data[start:end]
        else:
            return complement(self.data[start:end])

    def get_oriented(self, start, end, sense):
        if end < 0:
            return self.get_dummy(start, end, sense)
        elif start < 0:
            return self.get_dummy(start, 0, sense) + self.get_oriented(0, end, sense)
        if sense == "+":
            return self.data[start:end]
        else:
            try:
                return rev_comp(self.data[start:end])
            except KeyError:
                print start, end, sense

    def get_dummy(self, start, end, sense):
        return "N" * (end - start)

class Hit(object):
    def __init__(self):
        self.reads = []
        self.uniq = set()
        self.mapquals = []
        self.tissues = defaultdict(int)
        self.edits = []
        self.overlaps = []
        self.n_hits = []

    def add(self, read, A, B, dist, ov, n_hits, minmapscore, samples):
        self.reads.append(read)
        self.edits.append(dist)
        self.overlaps.append(ov)
        self.n_hits.append(n_hits)

        # Alignment Score - Secondbest hit score
        aopt = dict(A.tags)
        bopt = dict(B.tags)
        qA = aopt.get('AS') - aopt.get('XS', minmapscore)
        qB = bopt.get('AS') - bopt.get('XS', minmapscore)

        self.mapquals.append((qA + qB, qA, qB))

        for (prefix, tiss) in samples:
            if A.qname.startswith(prefix):
                self.tissues[tiss] += 1
                break

        self.uniq.add((read, tiss))
        self.uniq.add((rev_comp(read), tiss))

    def scores(self, chrom, start, end, sense, wiggle):
        n_reads = len(self.reads)
        n_uniq = len(self.uniq) / 2

        total_mq, best_qual_A, best_qual_B = sorted(self.mapquals, reverse=True)[0]

        wiggle = numpy.arange(-wiggle, wiggle + 1)

        spliced_at_begin = 0
        for x in wiggle:
            begin = (chrom, start + x, sense)
            if begin in loci:
                for l in loci[begin]:
                    spliced_at_begin += len(l.reads)

        spliced_at_end = 0
        for x in wiggle:
            ending = (chrom, end + x, sense)
            if ending in loci:
                for l in loci[ending]:
                    spliced_at_end += len(l.reads)

        tissues = sorted(self.tissues.keys())
        tiss_counts = [str(self.tissues[k]) for k in tissues]
        return (n_reads, n_uniq, best_qual_A, best_qual_B, spliced_at_begin, spliced_at_end, tissues, tiss_counts, min(self.edits), min(self.overlaps), min(self.n_hits))

def worker_init(genome_path):
    # Initialize a Track in each worker
    global WORKER_GENOME
    WORKER_GENOME = Track(genome_path, GenomeAccessor)

def worker_process(task):
    try:
        return worker_process_inner(task)
    except Exception:
        tb = traceback.format_exc()
        errfile = 'worker_error{}_{}.log'.format(os.getpid(), int(time.time()))
        with open(errfile, 'w') as f:
            f.write('task: {}\n\n'.format(repr(task)))
            f.write(tb)
        sys.stderr.write('Worker exception; logged to {}\n'.format(errfile))
        sys.stderr.flush()
        raise

def worker_process_inner(task):
    """
    task: (key, chrom, read, a_aend, b_pos, asize, margin, maxdist)
    Returns: (key, hits)
    """
    key, chrom, read, a_aend, b_pos, asize, margin, maxdist = task
    Seg = namedtuple('Seg', ['pos', 'aend'])
    A = Seg(pos=0, aend=a_aend)
    B = Seg(pos=b_pos, aend=0)
    hits = find_breakpoints(A, B, read, WORKER_GENOME, chrom, asize, margin, maxdist)
    return (key, hits)


def main():
    global N, loci, circs, splices  # make visible to output()

    usage = """

      bowtie2 anchors.qfa.gz | %prog > candidates.bed 2> candidates.reads

    """

    parser = OptionParser(usage=usage)
    parser.add_option("-S", "--system", dest="system", type=str, default="", help="model system database (alternatively use -G <path_to_genome_folder>)")
    parser.add_option("-G", "--genome", dest="genome", type=str, default="", help="path to genome. Point to folder with one fasta file for each chromosome.")
    parser.add_option("-a", "--anchor", dest="asize", type=int, default=20, help="anchor size (default=20)")
    parser.add_option("-w", "--wiggle", dest="wiggle", type=int, default=2, help="maximum nts a linear splice site may be away from circ to be counted as competitor (default=2)")
    parser.add_option("-m", "--margin", dest="margin", type=int, default=2, help="maximum nts the BP is allowed to reside inside an anchor (default=2)")
    parser.add_option("-d", "--maxdist", dest="maxdist", type=int, default=2, help="maximum mismatches (no indels) allowed in anchor extensions (default=2)")
    parser.add_option("-p", "--prefix", dest="prefix", default="", help="prefix to prepend to each name")
    parser.add_option("-q", "--min_uniq_qual", dest="min_uniq_qual", type=int, default=2, help="minimal uniqness (alignment score margin to next-best hit) for anchor alignments (default=2)")
    parser.add_option("-s", "--stats", dest="stats", default="runstats.log", help="where to put stats (default='runstats.log'")
    parser.add_option("-r", "--reads2samples", dest="reads2samples", default="", help="path to tab-separated two-column file with read-name prefix -> sample ID mapping")
    parser.add_option("-j", "--jobs", dest="jobs", type=int, default=mp.cpu_count(), help="number of worker processes (default=cpu_count())")
    parser.add_option("-b", "--batch-size", dest="batch_size", type=int, default=500, help="records per batch submitted to the pool (default=500)")

    options, args = parser.parse_args()

    # Genome source: prefer -G path (folder with chrom.fa files)
    if options.system:
        # Placeholder for system-based genome (if available in your env)
        # genome = getattr(sequence_data.systems, options.system).genome
        genome_path = options.genome  # fallback
    else:
        genome_path = options.genome

    if options.reads2samples:
        samples = [line.rstrip().split('\t') for line in file(options.reads2samples)]
    else:
        samples = []

    samples.append(('', 'unknown'))

    minmapscore = options.asize * (-2)

    loci = defaultdict(list)
    circs = defaultdict(Hit)
    splices = defaultdict(Hit)

    N = defaultdict(int)

    sam = pysam.Samfile('-', 'r')

    NUM_WORKERS = max(1, int(options.jobs))
    BATCH_SIZE = max(1, int(options.batch_size))

    pool = mp.Pool(processes=NUM_WORKERS, initializer=worker_init, initargs=(genome_path,))

    key_counter = 0

    try:
        # Process in batches to bound memory and improve throughput
        def process_batch(batch_tasks, pending_meta):
            # Map this batch; stream results as they complete
            for key, hits in pool.imap_unordered(worker_process, batch_tasks):
                orient, Aprime_seg, Bprime_seg, readp, chrom = pending_meta[key]
                if not hits:
                    if orient == 'circ':
                        N['circ_no_bp'] += 1
                    else:
                        N['splice_no_bp'] += 1
                    continue

                if orient == 'circ':
                    N['circ_reads'] += 1
                else:
                    N['spliced_reads'] += 1

                n_hits = len(hits)
                for h in hits:
                    dist, ov, chrom_h, start, end, sense = h
                    if orient == 'circ':
                        h_key = (chrom_h, start + 1, end - 1, sense)
                        circs[h_key].add(readp, Aprime_seg, Bprime_seg, dist, ov, n_hits, minmapscore, samples)
                    else:
                        h_key = (chrom_h, start, end, sense)
                        splices[h_key].add(readp, Aprime_seg, Bprime_seg, dist, ov, n_hits, minmapscore, samples)
                        loci[(chrom_h, start, sense)].append(splices[h_key])
                        loci[(chrom_h, end, sense)].append(splices[h_key])

        batch_tasks = []
        pending_meta = {}

        for A, B in grouper(2, sam):
            if A is None or B is None:
                continue
            N['total'] += 1
            if A.is_unmapped or B.is_unmapped:
                N['unmapped'] += 1
                continue
            if A.tid != B.tid:
                N['other_chrom'] += 1
                continue
            if A.is_reverse != B.is_reverse:
                N['other_strand'] += 1
                continue

            dist = B.pos - A.pos
            if numpy.abs(dist) < options.asize:
                N['overlapping_anchors'] += 1
                continue

            # Determine orientation branch (circ vs splice)
            if (A.is_reverse and dist > 0) or (not A.is_reverse and dist < 0):
                orient = 'circ'
            else:
                orient = 'splice'

            # Extract original read and chrom
            try:
                read = A.qname.split('__')[1]
            except Exception:
                read = ""
            chrom = sam.getrname(A.tid)

            # Prepare swapped A',B' for both worker compute and parent Hit.add semantics
            if A.is_reverse:
                # Swap and reverse complement read as in original code
                readp = rev_comp(read)
                Aprime_seg = B
                Bprime_seg = A
                a_aend = int(getattr(B, 'aend', B.pos))
                b_pos = int(A.pos)
            else:
                readp = read
                Aprime_seg = A
                Bprime_seg = B
                a_aend = int(getattr(A, 'aend', A.pos))
                b_pos = int(B.pos)

            # Build task and pending metadata
            key_counter += 1
            task = (key_counter, chrom, readp, a_aend, b_pos, int(options.asize), int(options.margin), int(options.maxdist))
            batch_tasks.append(task)
            pending_meta[key_counter] = (orient, Aprime_seg, Bprime_seg, readp, chrom)

            # Dispatch batch if full
            if len(batch_tasks) >= BATCH_SIZE:
                process_batch(batch_tasks, pending_meta)
                batch_tasks = []
                pending_meta = {}

        # Process remaining tasks
        if batch_tasks:
            process_batch(batch_tasks, pending_meta)

    except KeyboardInterrupt:
        pool.terminate()
        pool.join()
        sys.exit(1)
    finally:
        pool.close()
        pool.join()
        sam.close()

    # Write stats and outputs
    stats = file(options.stats, "w")
    stats.write(str(N) + "\n")

    # Note: output() expects (cand, prefix, suffix, wiggle, min_uniq_qual)
    output(circs, options.prefix, "circ", options.wiggle, options.min_uniq_qual)
    output(splices, options.prefix, "norm", options.wiggle, options.min_uniq_qual)

if __name__ == "__main__":
    main()
