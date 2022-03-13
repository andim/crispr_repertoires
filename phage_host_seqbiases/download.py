import os.path
import urllib.request
import pandas as pd

from pylib import *

datadir = 'data/'
path = datadir+'assembly_summary_refseq.txt' 
if not os.path.exists(path):
    url = r"ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
    print(url)
    urllib.request.urlretrieve(url, path)

refseq = pd.read_csv(datadir+'assembly_summary_refseq.txt', sep='\t', skiprows=1,
                     index_col=0, dtype={'taxid': str},
                     usecols=(0, 5, 7, 19))
refseq.index = [ind.split('.')[0] for ind in refseq.index]
for taxonid in ['562', '287', '28901']:
    path = datadir + taxonid + '.fna.gz'
    if not os.path.exists(path):
        try:
            ftpdir = refseq.loc['GCF_'+taxonid]['ftp_path']
            ftppath = ftpdir + '/' + ftpdir.split('/')[-1] + '_genomic.fna.gz'
            print('downloading', path)
            urllib.request.urlretrieve(ftppath, path)
        except KeyError:
            print('Failed: key not found, try using taxonid')
            seqs = refseq[refseq['taxid']==taxonid]
            if len(seqs) == 0:
                print('Failed again')
            for i, (accession, row) in enumerate(seqs.iterrows()):
                if len(seqs) > 1:
                    path = datadir + taxonid + '_' + str(i) + '.fna.gz'
                if not os.path.exists(path):
                    ftpdir = row['ftp_path']
                    ftppath = ftpdir + '/' + ftpdir.split('/')[-1] + '_genomic.fna.gz'
                    print('downloading', path)
                    urllib.request.urlretrieve(ftppath, path)
                if i > 1:
                    break

