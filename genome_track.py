import os, subprocess
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.io as pio

from mylib.gtf import GTF

gtf = GTF('../ref/genomic.gtf')

TE_PATH = '../Ribo-seq/TE'
RIBO_PATH = '../Ribo-seq/deseq2_protein_coding_first_tenth'
RNA_PATH = '../RNA-seq/deseq2_protein_coding_first_tenth'
IMG_PATH = '.'  # directory containing images (relative to project)
SAMPLES = 'BC'
CONTROL_GROUP='A'



RNA_COV='/fs/ess/PAS2967/dengyw144/rrnh/RNA-seq/Combined_norm_genome_cov'
RIBO_COV='/fs/ess/PAS2967/dengyw144/rrnh/Ribo-seq/coverage_norm_merged'


def generate_image(gene_id,output, sample):
    # print(gene_id)
    start, end, strand, chrom=gtf.get_gene_positions_strand_chrom(gene_id)
    ext=0
    start-=ext
    end+=ext
    ribo_cov_cg=os.path.join(RIBO_COV, f'Ribo-A{strand}.bedgraph')
    rna_cov_cg=os.path.join(RNA_COV, f'A-{'negative' if strand=='-' else 'positive'}-norm-coverage.bedgraph')
    
    ribo_cov_tg=os.path.join(RIBO_COV, f'Ribo-{sample}{strand}.bedgraph')
    rna_cov_tg=os.path.join(RNA_COV, f'{sample}-{'negative' if strand=='-' else 'positive'}-norm-coverage.bedgraph')
    
    cmd1=f'make_tracks_file --trackFiles {ribo_cov_cg} {rna_cov_cg} {ribo_cov_tg} {rna_cov_tg} genomic.bed -o tracks.ini'
    cmd2=f'pyGenomeTracks --tracks tracks.ini --region {chrom}:{start}-{end} --outFileName {output}/{gene_id}.png --dpi 300'
    r=subprocess.run([cmd1], shell=True, check=True,  capture_output=True, text=True)
    r=subprocess.run([cmd2], shell=True, check=True,  capture_output=True, text=True)


    
def plot_1_sample(rna_file, ribo_file, te_file, output, sample):
    os.makedirs(output, exist_ok=True)

    rna_df = pd.read_csv(rna_file, index_col=0)
    ribo_df = pd.read_csv(ribo_file, index_col=0)
    te_df = pd.read_csv(te_file, sep='\t', index_col=0)


    for g in te_df.index:
        if g not in rna_df.index or g not in ribo_df.index:
            continue

        x = rna_df.at[g, 'log2FoldChange']
        y = ribo_df.at[g, 'log2FoldChange']
        padj_x = rna_df.at[g, 'padj']
        padj_y = ribo_df.at[g, 'padj']
        padj_te = te_df.at[g, 'padj']

        if pd.isna(padj_x) or pd.isna(padj_y) or pd.isna(padj_te):
            continue
        generate_image(g, output, sample)
        

for s in SAMPLES:
    rna_file = os.path.join(RNA_PATH, f"{s}.csv")
    ribo_file = os.path.join(RIBO_PATH, f"{s}.csv")
    te_file = os.path.join(TE_PATH, f"result-{s}.tsv")
    output = f'coverage_{s}'
    plot_1_sample(rna_file, ribo_file, te_file, output, s)