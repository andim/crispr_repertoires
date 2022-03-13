import os
from collections import defaultdict
from functools import partial
from mimetypes import guess_type
import xml.etree.ElementTree as ET
import gzip
import re

import numpy as np
import scipy.special
import pandas as pd
from Bio import SeqIO

# Define path variables
repopath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
datadir = os.path.join(repopath, 'data/')

def parse_taxonids(xmlfile):
    tree = ET.parse(xmlfile)
    taxonids = []
    for taxon in tree.getroot()[0].findall('Taxon'):
        taxonids.append(taxon.find('TaxonId').text)
    return taxonids
 

def parse(xmlfile):
    tree = ET.parse(xmlfile)
    taxonids = []
    nloci = []
    genomelengths = []
    spacerlengthss = []
    for taxon in tree.getroot()[0].findall('Taxon'):
        taxonids.append(taxon.find('TaxonId').text)
        sequences = taxon.find('Sequences')
        nloci.append(int(sequences.find('SequenceCRISPRCount').text))
        length = 0
        spacerlengths = []
        for seq in sequences.findall('Sequence'):
            length += int(seq.find('GenomeLength').text)
            for crispr in seq.find('CRISPRs').findall('CRISPR'):
                for spacer in crispr.find('Spacers').findall('Spacer'):
                    spacerlengths.append(int(spacer.find('SpacerLength').text))
        spacerlengthss.append(spacerlengths)
        genomelengths.append(length)
    mean_spacerlengths = [np.mean(cls) if cls else np.nan for cls in spacerlengthss]
    median_spacerlengths = [np.median(cls) if cls else np.nan for cls in spacerlengthss]
    std_spacerlengths = [np.std(cls) if cls else np.nan for cls in spacerlengthss]
    nspacers = [len(spacerlengths) for spacerlengths in spacerlengthss]
    df = pd.DataFrame.from_dict(dict(taxonid=taxonids,
                                 nloci=nloci,
                                 genomelength=genomelengths,
                                 nspacer=nspacers,
                                 mean_spacerlength=mean_spacerlengths,
                                 median_spacerlength=median_spacerlengths,
                                 std_spacerlength=std_spacerlengths))
    return df, spacerlengthss


def fasta_iter(fasta_name, returnheader=True, returndescription=False):
    """
    Given a fasta file return a iterator over tuples of header, complete sequence.
    """
    if returnheader and returndescription:
        raise Exception('one of returnheader/returndescription needs to be False')
    if guess_type(fasta_name)[1] =='gzip':
        _open = partial(gzip.open, mode='rt')
    else:
        _open = open
    with _open(fasta_name) as f:
        fasta_sequences = SeqIO.parse(f, 'fasta')
        for fasta in fasta_sequences:
            if returndescription:
                yield fasta.description, str(fasta.seq)
            elif returnheader:
                yield fasta.id, str(fasta.seq)
            else:
                yield str(fasta.seq)


def entropy_grassberger(n, base=None):
    """"
    Estimate the entropy of a discrete distribution from counts per category

    n: array of counts 
    base: base in which to measure the entropy (default: nats)
    """
    N = np.sum(n)
    entropy = np.log(N) - np.sum(n*scipy.special.digamma(n+1e-20))/N
    if base:
        entropy /= np.log(base)
    return entropy


class Counter(defaultdict):

    def __init__(self, iterable, k, gap=0, **kwargs):
        """
        Counter class

        iterable: sequences or fasta filename
        k: int, kmer length
        gap: int, gap between first and subsequent letters
        """
        super(Counter, self).__init__(int)
        self.k = k
        self.gap = gap
        if isinstance(iterable, str):
            iterable = fasta_iter(iterable, returnheader=False)
        self.count(iterable, **kwargs)

    def count(self, iterable, **kwargs):
        for seq in iterable:
            count_kmers(seq, self.k, gap=self.gap, counter=self, **kwargs)

    def to_df(self, norm=True):
        """Convert a (kmer, count) dict to a pandas DataFrame
        """
        if norm:
            return pd.DataFrame(dict(seq=list(self.keys()), freq=normalize(self)))
        arr = np.array(list(self.values()), dtype=np.float)
        return pd.DataFrame(dict(seq=list(self.keys()), count=arr))

try:
    from clib import count_kmers
except ImportError:
    def count_kmers(string, k, counter=None, gap=0):
        """
        Count occurrence of kmers in a given string.
        """
        if counter is None:
            counter = defaultdict(int)
        for i in range(len(string)-k-gap+1):
            if gap:
                counter[string[i]+string[i+gap+1:i+k+gap]] += 1
            else:
                counter[string[i:i+k]] += 1
        return counter
        
def normalize(counter):
    "Given a (kmer, count) dict returns a normalized array of frequencies"
    arr = np.array(list(counter.values()), dtype=np.float)
    arr /= np.sum(arr)
    return arr


nucleotide_set = set(['A', 'C', 'G', 'T'])
def isvalidnt(string):
    "returns true if string is composed only of characters from the standard nucleotide alphabet"
    return all(c in nucleotide_set for c in string)

def to_kmers(seqs, k, return_index=False):
    """Generator yielding all possible kmers from a set of sequences.

    seqs: list of sequences
    k: length of kmer
    return_index: if true yield tuple seq, index
    """
    if return_index:
        for index, seq in enumerate(seqs):
            for i in range(len(seq)-k+1):
                s = seq[i:i+k]
                if isvalidnt(s):
                    yield s, index
    else:
        for seq in seqs:
            for i in range(len(seq)-k+1):
                s = seq[i:i+k]
                if isvalidnt(s):
                    yield s

def load_fasta_as_df(path):
    headers, seqs = list(zip(*[(h, seq) for h, seq in fasta_iter(path,
        returndescription=True, returnheader=False)]))
    df = pd.DataFrame(dict(header=headers, sequence=seqs))
    return df


