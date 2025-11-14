import os
from math import log
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from mylib.gtf import GTF
from mylib.feature import Feature
from mylib.utils import get_paths_ends_with_something, _auto_grouping

gtf=GTF('../reference/genomic.gtf', '../reference/GCF_000750555.1_ASM75055v1_genomic.fna')

def first_third(lst):
    return sorted(lst, reverse=True)[:len(lst)//10]


RNA=True
prefix='RNA' if RNA else 'Ribo'
DIR_FEATURE=f"../{prefix}-seq/feature"
FEATURE_PATHS=get_paths_ends_with_something(DIR_FEATURE, '.txt')
FEATURES=[Feature(p) for p in FEATURE_PATHS]
NAMES=[f.name for f in FEATURES]

groups= _auto_grouping(FEATURES, key=lambda x: str(x))

GENE_LENGTH={}
def distribution(samples):
    result={}
    for sample in samples:
        sample: Feature
        for gene in sample.gene_ids:
            # if not gtf.is_protein(gene):
            #     continue
            if gene not in GENE_LENGTH:
                GENE_LENGTH[gene]=gtf.id2length(gene)
            if gene not in result:
                result[gene]=[]
            result[gene].append(sample.gene_count[gene]/ GENE_LENGTH[gene])
    for g in result:
        result[g]=np.mean(result[g])*100
    return result


THRESHOLDS={}
GROUP_GENE_COUNT={} # group: {gene_id: count}

for group in groups:
    dist=distribution(groups[group])
    GROUP_GENE_COUNT[group]=dist
    to_plot=dist.values()
    THRESHOLDS[group]=  min(first_third(to_plot))

def good(gene_id, group):
    try:
        return GROUP_GENE_COUNT[group][gene_id] >= THRESHOLDS[group]
    except KeyError:
        return False
# breakpoint()
GOOD_GENES=[]
groups=[i for i in 'ABC']

for g in gtf.all_genes():
    if any(good(g, group) for group in groups):
        GOOD_GENES.append(g)

def filter(path, path_out):
    df = pd.read_csv(path, sep="\t")
    df = df[df['Entry'].isin(GOOD_GENES)]
    df.to_csv(path_out, sep="\t", index=False)

# filter('_TE/matrix-1.tsv', 'TE/matrix-1.tsv')
# filter('_TE/matrix-2.tsv', 'TE/matrix-2.tsv')

# df = pd.read_csv('your_file.csv', index_col=0)
# df.index.name = 'gene_id'

with open('10p.pkl', 'wb') as f:
    import pickle
    pickle.dump(GOOD_GENES, f)
breakpoint()