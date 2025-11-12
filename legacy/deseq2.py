import os
from itertools import combinations
import subprocess
import fcntl

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from mylib.gtf import GTF

gtf=GTF('../reference/genomic.gtf')
protein_coding_genes=gtf.all_genes()

DESEQ2_SCRIPT='/fs/ess/PAS2967/dengyw144/S21/deseq2.R'

RNA_PATH='../Ribo-seq/intermediate/feature'

DESEQ2_PATH='../Ribo-seq/deseq2_protein_coding'


def get_summary_files(folder_path):
    txt_files = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".summary"):
                txt_files.append(os.path.join(root, file))
    return txt_files

def get_txt_files(folder_path):
    txt_files = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".txt"):
                txt_files.append(os.path.join(root, file))
    return txt_files

class Feature:
    '''
        Class to handle 1 featureCount file
    '''
    name: str
    counts: list    # int
    gene_ids: list # string
    def __init__(self, path, ribo_seq=False):
        df = pd.read_csv(path, sep="\t", comment="#")
        FIELD_COUNTS=df.columns.tolist()[-1]
        self.counts=df[FIELD_COUNTS].tolist()
        self.gene_ids=df['Geneid'].tolist()
        self.name=path.split('/')[-1].split('.')[0]
        if ribo_seq:
            self.name='Ribo-'+self.name
    def scatter_plot(self, other, base_path):
        os.system(f'mkdir -p {base_path}')
        plt.scatter(self.counts, other.counts)
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig(os.path.join(base_path, f"{self.name}-{other.name}.png"), dpi=300)
        plt.cla()


rnas= get_txt_files(RNA_PATH)

features_RNA=[Feature(path) for path in rnas]


##### Create scatter plots
if 1:
    for f1, f2 in combinations(features_RNA,2):
        f1.scatter_plot(f2, '../Ribo-seq/scatter')
quit()

##### Prepare RiboDiff input
# CONTROL_GROUPS=['W1', 'W2']
# TREATED_GROUPS_1=['A'+str(i) for i in [1,2,3]]
# TREATED_GROUPS_2=['B'+str(i) for i in [1,2,3]]


def auto_grouping(features, control_group='A'):
    cg=[]
    tg={}
    for f in features:
        group=f.name[0]
        if group==control_group:
            cg.append(f.name)
        else:
            if group not in tg:
                tg[group]=[]
            tg[group].append(f.name)
    return tg, cg

print(auto_grouping(features_RNA))


# os.makedirs(DESEQ2_PATH, exist_ok=True)
def make_deseq2_input(i, tg, cg, features):
    tg=tg[i]
    result=np.array([features[0].gene_ids] + [f.counts for f in features if any(g in f.name for g in tg+cg)]).T

    df=pd.DataFrame(result)
    df=df[df[0].isin(protein_coding_genes)]
    df.to_csv(os.path.join(DESEQ2_PATH, f"matrix-{i}.csv"), sep=',' , index=False, header=['GeneID']+[f.name for f in features if any(g in f.name for g in tg+cg)])

    # Make meta
    meta=[]
    for f in features:
        if any(g in f.name for g in tg):
            condition='treated'
        elif any(g in f.name for g in cg):
            condition='control'
        else:
            continue
        

        meta.append([f.name, condition])
    
    df=pd.DataFrame(meta)
    df.to_csv(os.path.join(DESEQ2_PATH,f"meta-{i}.csv"), sep=',' , index=False, header=['id','dex'])




tg, cg=auto_grouping(features_RNA, control_group='A')
if 1:
    os.makedirs(DESEQ2_PATH, exist_ok=True)

    for i in tg:
        make_deseq2_input(i, tg, cg, features_RNA)
        command=f'Rscript {DESEQ2_SCRIPT} {DESEQ2_PATH} {i}'
        # breakpoint()
        r=subprocess.run([command], shell=True, check=True,  capture_output=True, text=True)
        with open('out.log', 'a') as f: 
            fcntl.flock(f, fcntl.LOCK_EX)
            f.write(command+'\n')
            if r.stdout is not None:
                f.write(r.stdout)
            if r.stderr is not None:
                f.write(r.stderr)
            f.write('\n')
            fcntl.flock(f, fcntl.LOCK_UN)


        

#############  Percentage of rrna in Ribo-seq
# RRNA_SUMMARY_FOLDER='Ribo-seq/n_rrna'
# os.makedirs(RRNA_SUMMARY_FOLDER, exist_ok=True)

# def summary2dict(path):
#     df = pd.read_csv(path, sep="\t",  index_col=0)
#     df = df.transpose()

#     # Convert the row to a dictionary (only one row after transposing)
#     data_dict = df.iloc[0].to_dict()
#     return data_dict

# def path2basename(path):
#     """
#     Extracts the base name of a file or directory from its path.
#     """
#     return os.path.basename(os.path.normpath(path)).split('.')[0]

# rrna_summaries=get_summary_files('Ribo-seq/feature_rrna')
# gene_summaries=get_summary_files('Ribo-seq/feature_renamed')

# rrna_dicts= {path2basename(path) :summary2dict(path) for path in rrna_summaries}
# gene_dicts= {path2basename(path) :summary2dict(path) for path in gene_summaries}
# results=[]
# for sample in rrna_dicts:
#     if sample not in gene_dicts:
#         raise Exception(f"Sample {sample} not found in gene summaries")
#     rrna=rrna_dicts[sample]
#     gene=gene_dicts[sample]
#     total_reads=gene['Assigned']
#     rrna_reads=rrna['Assigned']
#     unassigned=gene['Unassigned_Unmapped']+gene['Unassigned_NoFeatures']
#     perc_rrna=rrna_reads/total_reads
#     results.append([sample, total_reads, rrna_reads, unassigned, perc_rrna, unassigned/(total_reads+unassigned)])

# results= sorted(results, key=lambda x: x[0])  # Sort by Sample name
# df=pd.DataFrame(results, columns=['Sample', 'All_Assigned_Reads (gene+transcript)', 'rRNA, tRNA Reads', 'Unassigned_Reads', 'Percentage_rRNA, tRNA (among assigned)', 'Percentage_Unassigned'])
# df.to_csv(os.path.join(RRNA_SUMMARY_FOLDER, 'rrna_summary.csv'), index=False)
# print("rRNA summary saved to:", os.path.join(RRNA_SUMMARY_FOLDER, 'rrna_summary.csv'))