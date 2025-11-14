import os, subprocess, multiprocessing
import pandas as pd
import numpy as np
import ray

from mylib.gtf import GTF

gtf = GTF('reference/genomic.gtf')

N_CPU=multiprocessing.cpu_count()

SAMPLES = 'BC'
CONTROL_GROUP='A'



RNA_COV='/fs/ess/PAS2967/dengyw144/hbeak2/RNA-seq/bam'
RIBO_COV='/fs/ess/PAS2967/dengyw144/hbeak2/Ribo-seq/bam'
print('N_CPU: ', N_CPU)
ray.init(num_cpus=N_CPU)
os.makedirs('tmp', exist_ok=True)


GENE_IDS=gtf.all_genes()


def gene_id_to_args(gene_id):
    start, end, strand, chrom=gtf.get_gene_positions_strand_chrom(gene_id)
    ext=0
    start-=ext
    end+=ext
    return (start, end, strand, chrom)

GENE_ID_TO_ARGS={gene_id:gene_id_to_args(gene_id)  for gene_id in GENE_IDS}



@ray.remote(num_cpus=1)
def generate_image(gene_id,output, sample):
    if os.path.exists(f'{output}/{gene_id}.png'):
        return 0
    # print(gene_id)
    start, end, strand, chrom=GENE_ID_TO_ARGS[gene_id]

    ribo_cov_cg=os.path.join(RIBO_COV, f'Ribo-A.{"plus" if strand=='+' else 'minus'}.bedgraph')
    rna_cov_cg=os.path.join(RNA_COV, f'RNA-A.{"plus" if strand=='+' else 'minus'}.bedgraph')
    
    ribo_cov_tg=os.path.join(RIBO_COV, f'Ribo-{sample}.{"plus" if strand=='+' else 'minus'}.bedgraph')
    rna_cov_tg=os.path.join(RNA_COV, f'RNA-{sample}.{"plus" if strand=='+' else 'minus'}.bedgraph')
    
    cmd1=f'make_tracks_file --trackFiles {ribo_cov_cg} {rna_cov_cg} {ribo_cov_tg} {rna_cov_tg} -o tmp/{output}.{gene_id}.ini'
    cmd2=f'pyGenomeTracks --tracks tmp/{output}.{gene_id}.ini --region {chrom}:{start}-{end} --outFileName {output}/{gene_id}.png --dpi 300'
    r=subprocess.run([cmd1], shell=True, check=True,  capture_output=True, text=True)
    r=subprocess.run([cmd2], shell=True, check=True,  capture_output=True, text=True)
    return 0


    

        
futures_result=[]

for s in SAMPLES:
    output = f'coverage_{CONTROL_GROUP}{s}'
    os.makedirs(output, exist_ok=True)
    for gene_id in GENE_IDS:
        r=generate_image.remote(gene_id, output, s)
        futures_result.append(r)

ray.get(futures_result)

ray.shutdown()