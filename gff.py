#!/usr/bin/env python
from optparse import OptionParser
import os, sys

################################################################################
# gff
#
# Methods for working with Gene Feature Format files.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <cmd> <gff file>'
    parser = OptionParser(usage)
    parser.add_option('-d', dest='downstream', type='int', default=200, help='Downstream bp for promoters [Default: %default]')
    parser.add_option('-o', dest='out_file')
    parser.add_option('-u', dest='upstream', type='int', default=2000, help='Upstream bp for promoters [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide command and gff file')
    else:
        cmd = args[0]
        gff_file = args[1]

    if cmd.lower() in ['utr','utrs']:
        utrs(gff_file, options.out_file)
    elif cmd.lower() in ['ss','splice']:
        splice_sites(gff_file, output_file=options.out_file)
    elif cmd.lower() in ['intron','introns']:
        introns(gff_file, output_file=options.out_file)
    elif cmd.lower() in ['promoter','promoters']:
        promoters(gff_file, options.upstream, options.downstream, output_file=options.out_file)
    elif cmd.lower() in ['span','spans']:
        span_gene(gff_file, output_file=options.out_file)


################################################################################
# extend
#
# Given a gtf file, extend the genes on each side.
################################################################################
def extend(gtf_file, extend_5p=2000, extend_3p=2000, output_file=None):
    if not output_file:
        gtf_base = os.path.splitext(gtf_file)[0]
        output_file = '%s_ext.gtf' % gtf_base
    out = open(output_file, 'w')

    genes = read_genes(gtf_file)
    source = open(gtf_file).readline().split()[1]

    for gid in genes:
        g = genes[gid]

        if g.strand == '+':
            # print 5' extension
            if extend_5p > 0 and g.exons[0].start-extend_5p >= 1:
                cols = [g.chrom, source, '5p', str(g.exons[0].start-extend_5p), str(g.exons[0].start-1), '.', g.strand, '.', kv_gtf(g.kv)]
                print >> out, '\t'.join(cols)

            # print exons
            for ex in g.exons:
                cols = [g.chrom, source, 'exon', str(ex.start), str(ex.end), '.', g.strand, '.', kv_gtf(g.kv)]
                print >> out, '\t'.join(cols)

            # print 3' extension
            if extend_3p > 0:
                cols = [g.chrom, source, '3p', str(g.exons[-1].end+1), str(g.exons[-1].end+extend_3p), '.', g.strand, '.', kv_gtf(g.kv)]
                print >> out, '\t'.join(cols)

        else:
            # print 3' extension
            if extend_3p > 0 and g.exons[0].start-extend_3p >= 1:
                cols = [g.chrom, source, '3p', str(g.exons[0].start-extend_3p), str(g.exons[0].start-1), '.', g.strand, '.', kv_gtf(g.kv)]
                print >> out, '\t'.join(cols)

            # print exons
            for ex in g.exons:
                cols = [g.chrom, source, 'exon', str(ex.start), str(ex.end), '.', g.strand, '.', kv_gtf(g.kv)]
                print >> out, '\t'.join(cols)

            # print 5' extension
            if extend_5p > 0:
                cols = [g.chrom, source, '5p', str(g.exons[-1].end+1), str(g.exons[-1].end+extend_5p), '.', g.strand, '.', kv_gtf(g.kv)]
                print >> out, '\t'.join(cols)

    out.close()

################################################################################
# g2t
#
# Given a gtf file, return a mapping of gene_id's to sets of transcript_id's
################################################################################
def g2t(gtf_file='%s/research/common/data/lncrna/lnc_catalog.gtf' % os.environ['HOME']):
    d = {}
    for line in open(gtf_file):
        a = line.split('\t')
        kv = gtf_kv(a[8])
        d.setdefault(kv['gene_id'],set()).add(kv['transcript_id'])
    return d


################################################################################
# gtf_kv
#
# Convert the last gtf section of key/value pairs into a dict.
################################################################################
def gtf_kv(s):
    d = {}

    a = s.split(';')
    for key_val in a:
        if key_val.strip():
            eq_i = key_val.find('=')
            if eq_i != -1 and key_val[eq_i-1] != '"':
                kvs = key_val.split('=')
            else:
                kvs = key_val.split()

            if len(kvs) == 2:
                key = kvs[0]
                if kvs[1][0] == '"' and kvs[1][-1] == '"':
                    val = kvs[1].strip()[1:-1]
                else:
                    val = kvs[1].strip()
                d[key] = val

    return d


################################################################################
# introns
#
# Given a gtf file, find the introns
################################################################################
def introns(gtf_file, output_file=None):
    if not output_file:
        gtf_base = os.path.splitext(gtf_file)[0]
        output_file = '%s_introns.gtf' % gtf_base
    out = open(output_file, 'w')

    genes = read_genes(gtf_file)
    source = open(gtf_file).readline().split()[1]

    for gid in genes:
        g = genes[gid]

        for i in range(1,len(g.exons)):
            cols = [g.chrom, source, 'intron', str(g.exons[i-1].end+1), str(g.exons[i].start-1), '.', g.strand, '.', kv_gtf(g.kv)]
            print >> out, '\t'.join(cols)

    out.close()


################################################################################
# kv_gtf
#
# Convert a kv hash to str gtf representation.
################################################################################
def kv_gtf(d):
    s = ''
    for key in sorted(d.keys()):
        s += '%s "%s"; ' % (key,d[key])
    return s


################################################################################
# promoters
#
# Given a gtf file, find the promoters
################################################################################
def promoters(gtf_file, promoter_up=2000, promoter_down=0, output_file=None, guess_strand=False):
    if not output_file:
        gtf_base = os.path.splitext(gtf_file)[0]
        output_file = '%s_prom.gtf' % gtf_base
    out = open(output_file, 'w')

    genes = read_genes(gtf_file)
    source = open(gtf_file).readline().split()[1]

    for gid in genes:
        g = genes[gid]

        if g.strand == '+':
            if g.exons[0].start-promoter_up >= 1:
                cols = [g.chrom, source, 'promoter', str(g.exons[0].start-promoter_up), str(g.exons[0].start+promoter_down), '.', g.strand, '.', kv_gtf(g.kv)]
            else:
                print >> sys.stderr, 'WARNING: %s discluded for nearness to chromosome end' % gid
                cols = []
        elif g.strand == '-':
            if g.exons[-1].end-promoter_down >= 1:
                cols = [g.chrom, source, 'promoter', str(g.exons[-1].end-promoter_down), str(g.exons[-1].end+promoter_up), '.', g.strand, '.', kv_gtf(g.kv)]
            else:
                print >> sys.stderr, 'WARNING: %s discluded for nearness to chromosome end' % gid
                cols = []
        else:
            if guess_strand:
                print >> sys.stderr, 'WARNING: %s guessing forward strand' % gid
                cols = [g.chrom, source, 'promoter', str(g.exons[0].start-promoter_up), str(g.exons[0].start+promoter_down), '.', '+', '.', kv_gtf(g.kv)]
            else:
                print >> sys.stderr, 'WARNING: %s discluded for lack of strand' % gid
                cols = []

        if cols:
            print >> out, '\t'.join(cols)

    out.close()


################################################################################
# span_gene
#
# Given a gtf file, make a new gtf file with a single entry per gene covering
# the entire span of the gene.
################################################################################
def span_gene(gtf_file, output_file=None):
    if not output_file:
        gtf_base = os.path.splitext(gtf_file)[0]
        output_file = '%s_span.gtf' % gtf_base

    # get isoforms
    transcripts = read_genes(gtf_file, key_id='transcript_id', sort=True)

    # store gene info
    gene_regions = {}
    for tid in transcripts:
        tx = transcripts[tid]
        gid = tx.kv['gene_id']
        if gid not in gene_regions:
            gene_regions[gid] = [tx.chrom, tx.exons[0].start, tx.exons[-1].end, tx.strand]
        else:
            gene_regions[gid][1] = min(gene_regions[gid][1], tx.exons[0].start)
            gene_regions[gid][2] = max(gene_regions[gid][2], tx.exons[-1].end)

    # print
    source = open(gtf_file).readline().split()[1]
    out = open(output_file, 'w')
    for gid in gene_regions:
        g = gene_regions[gid]
        cols = [g[0], source, 'span', str(g[1]), str(g[2]), '.', g[3], '.', kv_gtf({'gene_id':gid})]
        print >> out, '\t'.join(cols)
    out.close()


################################################################################
# span_transcript
#
# Given a gtf file, make a new gtf file with a single entry per transcript
# covering the entire span of the transcript.
################################################################################
def span_transcript(gtf_file, output_file=None):
    if not output_file:
        gtf_base = os.path.splitext(gtf_file)[0]
        output_file = '%s_span.gtf' % gtf_base

    # get isoforms
    transcripts = read_genes(gtf_file, key_id='transcript_id', sort=True)

    # print span
    source = open(gtf_file).readline().split()[1]
    out = open(output_file, 'w')
    for tid in transcripts:
        tx = transcripts[tid]
        cols = [tx.chrom, source, 'span', str(tx.exons[0].start), str(tx.exons[-1].end), '.', tx.strand, '.', kv_gtf(tx.kv)]
        print >> out, '\t'.join(cols)
    out.close()


################################################################################
# splice_sites
#
# Given a gtf file, find the splice sites
################################################################################
def splice_sites(gtf_file, exon=0, intron=2, output_file=None):
    if not output_file:
        gtf_base = os.path.splitext(gtf_file)[0]
        output_file = '%s_splice.gff' % gtf_base
    out = open(output_file, 'w')

    genes = read_genes(gtf_file)
    source = open(gtf_file).readline().split()[1]

    for gid in genes:
        g = genes[gid]

        for i in range(len(g.exons)-1):
            if g.strand == '+':
                # donor
                cols = [g.chrom, source, 'donor', str(g.exons[i].end-exon+1), str(g.exons[i].end+intron), '.', g.strand, '.', kv_gtf(g.kv)]
                print >> out, '\t'.join(cols)

                # acceptor
                cols = [g.chrom, source, 'acceptor', str(g.exons[i+1].start-intron), str(g.exons[i+1].start+exon-1), '.', g.strand, '.', kv_gtf(g.kv)]
                print >> out, '\t'.join(cols)
            else:
               # donor
               cols = [g.chrom, source, 'donor', str(g.exons[i+1].start-intron), str(g.exons[i+1].start+exon-1), '.', g.strand, '.', kv_gtf(g.kv)]
               print >> out, '\t'.join(cols)
           
               # acceptor
               cols = [g.chrom, source, 'acceptor', str(g.exons[i].end-exon+1), str(g.exons[i].end+intron), '.', g.strand, '.', kv_gtf(g.kv)]
               print >> out, '\t'.join(cols)

    out.close()


################################################################################
# read_genes
#
# Parse a gtf file and return a set of Gene objects in a hash keyed by the
# id given.
#
# Note: assumes exons only.
################################################################################
def read_genes(gtf_file, key_id='transcript_id', sort=True):
    genes = {}
    for line in open(gtf_file):
        a = line.split('\t')

        if a[2] in ['exon','CDS']:
            kv = gtf_kv(a[8])
            if not kv[key_id] in genes:
                genes[kv[key_id]] = Gene(a[0], a[6], kv)

            if a[2] == 'exon':
                genes[kv[key_id]].add_exon(int(a[3]), int(a[4]), sort=sort)
            elif a[2] == 'CDS':
                genes[kv[key_id]].add_cds(int(a[3]), int(a[4]), sort=sort)

    return genes


################################################################################
# three_prime
#
# Given a gtf file, return a section surrounding the 3' end.
################################################################################
def three_prime(gtf_file, upstream=0, downstream=2000, output_file=None):
    if not output_file:
        gtf_base = os.path.splitext(gtf_file)[0]
        output_file = '%s_3p.gtf' % gtf_base
    out = open(output_file, 'w')

    genes = read_genes(gtf_file)
    source = open(gtf_file).readline().split()[1]

    for gid in genes:
        g = genes[gid]

        if g.strand == '+' and g.exons[-1].end-upstream >= 1:
            cols = [g.chrom, source, '3p', str(g.exons[-1].end-upstream), str(g.exons[-1].end+downstream), '.', g.strand, '.', kv_gtf(g.kv)]
        elif g.strand == '-' and g.exons[0].start-downstream >= 1:
            cols = [g.chrom, source, '3p', str(g.exons[0].start-downstream), str(g.exons[0].start+upstream), '.', g.strand, '.', kv_gtf(g.kv)]
        else:
            cols = []
        
        if cols:
            print >> out, '\t'.join(cols)

    out.close()


################################################################################
# t2g
#
# Given a gtf file, return a mapping of transcript to gene id's
################################################################################
def t2g(gtf_file='%s/research/common/data/lncrna/lnc_catalog.gtf' % os.environ['HOME']):
    d = {}
    for line in open(gtf_file):
        a = line.split('\t')
        kv = gtf_kv(a[8])
        d[kv['transcript_id']] = kv['gene_id']
    return d


################################################################################
# utrs
#
# Take a gtf file with exons and CDS annotated and return a gtf of the UTRs
################################################################################
def utrs(gtf_file, output_file=None):
    if not output_file:
        gtf_base = os.path.splitext(gtf_file)[0]
        output_file = '%s_utrs.gtf' % gtf_base
    out = open(output_file, 'w')

    genes = read_genes(gtf_file)
    source = open(gtf_file).readline().split()[1]

    for gid in genes:
        g = genes[gid]
        
        # match up exons and CDS
        c = 0
        for e in range(len(g.exons)):
            # left initial
            if g.exons[e].end < g.cds[c].start:
                utr_label = (g.strand == '+')*'5\'UTR' + (g.strand == '-')*'3\'UTR'
                cols = [g.chrom, source, utr_label, str(g.exons[e].start), str(g.exons[e].end), '.', g.strand, '.', kv_gtf(g.kv)]
                print >> out, '\t'.join(cols)

            # right initial
            elif g.cds[c].end < g.exons[e].start:
                utr_label = (g.strand == '+')*'3\'UTR' + (g.strand == '-')*'5\'UTR'
                cols = [g.chrom, source, utr_label, str(g.exons[e].start), str(g.exons[e].end), '.', g.strand, '.', kv_gtf(g.kv)]
                print >> out, '\t'.join(cols)

            # overlap
            else:
                # left overlap
                if g.exons[e].start < g.cds[c].start:
                    utr_label = (g.strand == '+')*'5\'UTR' + (g.strand == '-')*'3\'UTR'
                    cols = [g.chrom, source, utr_label, str(g.exons[e].start), str(g.cds[c].start-1), '.', g.strand, '.', kv_gtf(g.kv)]
                    print >> out, '\t'.join(cols)

                # right overlap
                if g.cds[c].end < g.exons[e].end:
                    utr_label = (g.strand == '+')*'3\'UTR' + (g.strand == '-')*'5\'UTR'
                    cols = [g.chrom, source, utr_label, str(g.cds[c].end+1), str(g.exons[e].end), '.', g.strand, '.', kv_gtf(g.kv)]
                    print >> out, '\t'.join(cols)

                c = min(c+1,len(g.cds)-1)

    out.close()


################################################################################
# Gene
################################################################################
class Gene:
    def __init__(self, chrom, strand, kv):
        self.chrom = chrom
        self.strand = strand
        self.kv = kv
        self.exons = []
        self.cds = []

    def add_cds(self, start, end, sort=True):
        self.cds.append(Exon(start,end))
        if sort and len(self.cds) > 1 and self.cds[-2].end > start:
            #print >> sys.stderr, 'CDS are not sorted - %s %d %d' % (self.chrom,self.start,self.end)
            self.cds.sort()

    def add_exon(self, start, end, sort=True):
        self.exons.append(Exon(start,end))
        if sort and len(self.exons) > 1 and self.exons[-2].end > start:
            #print >> sys.stderr, 'Warning: exons are not sorted - %s' % kv_gtf(self.kv)
            self.exons.sort()

    def __str__(self):
        return '%s %s %s %s' % (self.chrom, self.strand, kv_gtf(self.kv), ','.join([ex.__str__() for ex in self.exons]))


################################################################################
# Exon
################################################################################
class Exon:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __cmp__(self, x):
        if self.start < x.start:
            return -1
        elif self.start > x.start:
            return 1
        else:
            return 0

    def __str__(self):
        return 'exon(%d-%d)' % (self.start,self.end)
        

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
