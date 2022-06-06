#!/usr/bin/env python

import argparse
import os
import re
import sys
import gzip
import tempfile
import time
from math import exp, ldexp, log, ceil
from mmh3 import hash64
from random import randrange
from statistics import mean, stdev
from subprocess import call

__author__ = "Ka Ming Nip"
__copyright__ = "Copyright 2018, Canada's Michael Smith Genome Sciences Centre"


def get_num_kmers(seq, k):
    #calculates the number of kmers
    return len(seq) - k + 1
# enddef

def kmers_stream(seq, k):
    #iterate kmers of the sequence
    for i in range(get_num_kmers(seq, k)):
        yield seq[i:i + k]
    # endfor
# enddef

def kmers_set(seq, k):
    #return the unique set of kmers
    return set(kmers_stream(seq, k=k))
# enddef

def kmers_stream_multi(seqs, k):
    for s in seqs:
        for kmer in kmers_stream(s, k):
            yield kmer
        #endfor
    #endfor
#enddef

def kmers_set_multi(seqs, k):
    #return the unique set of kmers
    return set(kmers_stream_multi(seqs, k=k))
# enddef

complements_trans = str.maketrans('ACGTacgt', 'TGCATGCA')

def reverse_complement(seq):
    #create reverse complement
    return seq[::-1].translate(complements_trans)
# enddef

def non_stranded_hash(string):
    return min(hash64(string, signed=False), hash64(reverse_complement(string), signed=False))
# enddef

def fasta_seq_stream(fh):
    #read through fasta file
    name = None
    seq = []

    for line in fh:
        firstchar = line[0]

        if firstchar == '>':
            # the previous seq
            if len(seq) > 0:
                yield (name, ''.join(seq))
            # endif

            seq = []
            name = line.rstrip().split(maxsplit=1)[0][1:]
        else:
            seq.append(line.rstrip())
        # enddef
    # enddef

    # the last seq of the file
    yield (name, ''.join(seq))

    fh.close()
# enddef

def make_acgt_re_obj(min_len=1):
    #regular expression for ACGT string of the specified min length
    return re.compile('[ACTG]{%d,}' % min_len, re.IGNORECASE)
# enddef

phred33 = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJ'

def fastq_seq_stream(fh, q_thres=3, min_len=1):
    qual_re_obj = re.compile('[%s]{%d,}' % (phred33[q_thres:], min_len))

    seq_re_obj = make_acgt_re_obj(min_len)

    try:
        while 1:
            assert next(fh)[0] == '@'  # name
            seq = next(fh)  # sequence
            assert next(fh)[0] == '+'  # comment

            for m in qual_re_obj.finditer(next(fh)):
                for s in seq_re_obj.findall(seq[m.start():m.end()]):
                    yield s.upper()
                # endfor
            # endfor
        # endwhile
    except StopIteration:
        pass
    # endtry

    fh.close()
# enddef

def get_cardinality(k, outprefix, read_files, verbose=True):
    # run ntCard, eg. ntcard -k 50 -p prefix reads.fq
    cmd = 'ntcard -k %d -p %s %s' % (k, outprefix, ' '.join(read_files))
    
    if verbose:
        sys.stderr.write('Running command:\n%s\n...\n' % cmd)
    
    if call(cmd, shell=True) != 0:
        sys.stderr.write('ERROR: ntCard encounter an error!\n')
        sys.exit(1)
    #endif

    # histogram file path, eg. prefix_k50.hist
    hist = outprefix + '_k' + str(k) + '.hist'
    
    if verbose:
        sys.stderr.write('Parsing histogram `%s`...\n' % hist)
    
    # extract number of unique kmers from histogram file
    num_kmers = 0
    with open(hist) as fh:
        for line in fh:
            if len(line) > 0:
                items = line.strip().split('\t')
                if items[0] == 'F0':
                    num_kmers = int(items[1])
                    break
                #endif
            #endif
        #endfor
    #endwith
    
    if verbose:
        sys.stderr.write('%d unique k-mers found\n' % num_kmers)
    
    # remove histogram file
    os.remove(hist)
    
    return num_kmers
#enddef

def get_stats(sorted_list):
    if sorted_list is None:
        return None
    # endif

    l = len(sorted_list)

    if l == 0:
        return None
    # endif

    minimum = sorted_list[0]
    maximum = sorted_list[-1]

    q1_i = l // 4
    m_i = l // 2
    q3_i = l // 2 + l // 4

    median = 0
    if l % 2:
        median = (sorted_list[m_i - 1] + sorted_list[m_i]) // 2
    else:
        median = sorted_list[m_i]
    # endif

    q1 = 0
    q3 = 0
    if l % 4:
        q1 = (sorted_list[q1_i - 1] + sorted_list[q1_i]) // 2
        q3 = (sorted_list[q3_i - 1] + sorted_list[q3_i]) // 2
    else:
        q1 = sorted_list[q1]
        q3 = sorted_list[q3]
    # endif

    return [minimum, q1, median, q3, maximum]
# enddef

class CountingBloomFilter:
    """A probablistic logarithmic counting bloom filter
    """

    # 5 exponent bits and 3 mantissa bits in minifloat represented as uint8
    mantissa = 3
    manti_mask = 0xFF >> (8 - mantissa)
    add_mask = 0x80 >> (7 - mantissa)
    name = 'CountingBloomFilter'

    def __init__(self, k, size=0, hash_count=2):
        self.size = size
        self.counts_array = bytearray(b'\x00' * size)
        self.hash_count = hash_count
        self.k = k
    # enddef

    def serialize(self, path):
        with open(path + '.desc', 'wt') as f:
            f.write('type:%s\n' % self.name)
            f.write('size:%d\n' % self.size)
            f.write('hash_count:%d\n' % self.hash_count)
            f.write('k:%d\n' % self.k)
            f.write('fdr:%s\n' % self.get_fdr())
        # endwith

        with open(path, 'wb') as f:
            f.write(self.counts_array)
        # endwith
    # enddef

    def deserialize(self, path):
        with open(path + '.desc', 'rt') as f:
            for line in f:
                key, val = line.rstrip().split(':')
                if key == 'type':
                    assert val == self.name
                elif key == 'size':
                    self.size = int(val)
                elif key == 'hash_count':
                    self.hash_count = int(val)
                elif key == 'k':
                    self.k = int(val)
                # endif
            # endfor
        # endwith

        self.counts_array = bytearray(b'\x00' * self.size)

        with open(path, 'rb') as f:
            f.readinto(self.counts_array)
        # endwith

        assert self.size == len(self.counts_array)
    # enddef

    def get_index(self, val):
        return val % self.size
    # enddef

    def get_fpr(self):
        # (n/m) ^ h
        h = self.hash_count
        m = self.size
        n = m - self.counts_array.count(b'\x00')
        
        return (n / m) ** h
    # enddef

    def add(self, string, stranded=False):
        if stranded:
            self.add_hash(hash64(string, signed=False))
        else:
            self.add_hash(non_stranded_hash(string))
        # endif
    # enddef

    def add_hash(self, hash_vals):
        # find min val
        current_min = self.counts_array[self.get_index(hash_vals[0])]
        for i in range(1, self.hash_count):
            c = self.counts_array[self.get_index(hash_vals[i])]

            if c < current_min:
                current_min = c
            # endif
            if current_min == 0:
                break
            # endif
        # endfor

        # increment count
        if (current_min <= self.manti_mask
                or (current_min < 255
                    and randrange(0xFFFFFFFF) % (1 << ((current_min >> self.mantissa) - 1)) == 0)):
            updated_min = current_min + 1

            for i in range(self.hash_count):
                if self.counts_array[self.get_index(hash_vals[i])] == current_min:
                    self.counts_array[self.get_index(hash_vals[i])] = updated_min
                # endif
            # endfor
        # endif
    # enddef

    def add_stream(self, stream, stranded=False):
        try:
            if stranded:
                while 1:
                    self.add_hash(hash64(next(stream), signed=False))
                # endwhile
            else:
                while 1:
                    self.add_hash(non_stranded_hash(next(stream)))
                # endwhile
            # endif
        except StopIteration:
            pass
        # endtry
    # enddef

    def lookup_hash(self, hash_vals):
        # find min val
        current_min = self.counts_array[self.get_index(hash_vals[0])]
        for i in range(1, self.hash_count):
            c = self.counts_array[self.get_index(hash_vals[i])]
            if c < current_min:
                current_min = c
            # endif
            if current_min == 0:
                return 0
            # endif
        # endfor

        if current_min <= self.manti_mask:
            return float(current_min)
        # endif

        return float(ldexp((current_min & self.manti_mask) | self.add_mask, (current_min >> self.mantissa) - 1))
    # enddef

    def lookup(self, string, stranded=False):
        hash_vals = [0] * self.hash_count

        if stranded:
            return self.lookup_hash(hash64(string, signed=False))
        else:
            return self.lookup_hash(non_stranded_hash(string))
        # endif
    # enddef

    def lookup_stream(self, stream, stranded=False):
        if stranded:
            for string in stream:
                yield self.lookup_hash(hash64(string, signed=False))
            # endfor
        else:
            for string in stream:
                yield self.lookup_hash(non_stranded_hash(string))
            # endfor
        # endif
    # enddef

    def lookup_seq(self, seq, acgt_re_obj, stranded=False, unique_kmers=False):
        kmerize = kmers_set_multi if unique_kmers else kmers_stream_multi
        for c in self.lookup_stream(kmerize(acgt_re_obj.findall(seq.upper()), k=k), stranded=stranded):
            yield c
        # endfor
    # enddef
# endclass

parser = argparse.ArgumentParser(description='Returns the distribution of k-mer counts of input sequences')
parser.add_argument('-s', '--seq', dest='sequences', metavar='FASTA', type=str, nargs='+',
                    help='sequence FASTA files')
parser.add_argument('-r', '--reads', dest='reads', metavar='FASTQ', type=str, nargs='+',
                    help='read FASTQ files')
parser.add_argument('-k', '--kmer', dest='k', metavar='INT', type=int, default=31,
                    help='kmer size [%(default)s]')
#parser.add_argument('-w', '--window', dest='window', metavar='INT', type=int, default=5,
#                    help='number of kmers in sliding window [%(default)s]')
parser.add_argument('--hash', dest='hash_count', metavar='INT', type=int, default=2,
                    help='number of hash values per kmer [%(default)s]')
parser.add_argument('--fpr', dest='fpr', metavar='FLOAT', type=float, default=0.05,
                    help='max false positive rate allowed [%(default)s]')
parser.add_argument('--ntcard', dest='ntcard', action='store_true',
                    help='use ntCard to help configure the Bloom filter size')
parser.add_argument('--noheader', dest='no_header', action='store_true', default=False,
                    help='do not print column header')
args = parser.parse_args()

#########################
# Execute the script here
#########################

for p in args.sequences:
    if not os.path.isfile(p):
        sys.stderr.write('ERROR: Input sequence file does not exist `%s`\n' % p)
        sys.exit(1)
    # endif
# endfor

files_num_bytes = 0
for p in args.reads:
    if not os.path.isfile(p):
        sys.stderr.write('ERROR: Input read file does not exist `%s`\n' % p)
        sys.exit(1)
    else:
        files_num_bytes += os.path.getsize(p)
    # endif
# endfor

hash_count = args.hash_count
k = args.k
#w = args.window
max_fpr = args.fpr
cbf_num_bytes = files_num_bytes

if args.ntcard:
    outprefix = os.path.join(tempfile.gettempdir(), 
                             tempfile.gettempprefix()) + '-' + os.path.basename(args.reads[0] + '-' + str(time.time()))
    
    num_unique_kmers = get_cardinality(k, outprefix, args.reads)
    
    # calculate Bloom filter size
    r = (-hash_count) / log(1 - exp(log(max_fpr) / hash_count));
    cbf_num_bytes = ceil(num_unique_kmers * r);
#endif

cbf = CountingBloomFilter(k=k, size=cbf_num_bytes, hash_count=hash_count)

while True:
    for p in args.reads:
        sys.stderr.write("Parsing `%s` ...\n" % p)
        if p.endswith('.gz'):
            fh = gzip.open(p, 'rt')
        else:
            fh = open(p, 'rt')
        # endif

        for seq in fastq_seq_stream(fh, min_len=k):
            cbf.add_stream(kmers_stream(seq, k=k))
        # endfor

        fh.close()
    # endfor

    fpr = cbf.get_fpr()

    if args.ntcard or fpr <= max_fpr:
        break
    #endif

    sys.stderr.write("FPR is higher than threshold (ie. %f > %f), resizing Bloom filter...\n" % (fpr, max_fpr))

    # FPR is too high, double the Bloom filter size
    cbf_num_bytes *= 2
    cbf = CountingBloomFilter(k=k, size=cbf_num_bytes, hash_count=hash_count)
#endwhile

sys.stderr.write("Bloom filter size: %f GB\n" % (cbf_num_bytes / float(1024 ** 3)))
sys.stderr.write("FPR: %f\n" % fpr)

acgt_re_obj = make_acgt_re_obj(k)

if not args.no_header:
    print("#name min q1 median q3 max mean stdev")
# endif

for p in args.sequences:
    sys.stderr.write("Parsing `%s` ...\n" % p)
    if p.endswith('.gz'):
        fh = gzip.open(p, 'rt')
    else:
        fh = open(p, 'rt')
    # endif

    for name, seq in fasta_seq_stream(fh):
        #counts = []
        slider = []
        #counts2 = []
        slider_uniq = []
        #print(cbf.lookup_seq(seq, acgt_re_obj))

        for x in cbf.lookup_seq(seq, acgt_re_obj):
            slider.append(x)
        #print(slider)

        for x in cbf.lookup_seq(seq, acgt_re_obj, unique_kmers=True):
            slider_uniq.append(x)
            #if len(slider) >= w:
            #    counts.append(mean(slider))
            #    slider.pop(0)
            #endif
        # endfor
        slider.sort()
        slider_uniq.sort()

        print(name, ' '.join([str(x) for x in get_stats(slider)]), mean(slider), stdev(slider))
        print(name + "_unique", ' '.join([str(x) for x in get_stats(slider_uniq)]), mean(slider_uniq), stdev(slider_uniq))

    # endfor
# endfor
