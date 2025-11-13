import os
import multiprocessing
import subprocess
import fcntl
import ray
import pysam
from mylib.utils import get_paths_ends_with_something
from func import *

BASE_DIR='Ribo-seq'
# BASE_DIR='test_data'

if __name__=='__main__':
    print('N_CPU: ', N_CPU)
    ray.init(num_cpus=N_CPU)
    
    # retrieve sample names
    SAMPLES=[i+j for i in 'ABC' for j in '123']
    
    DEDUP_DIR=os.path.join(BASE_DIR, 'dedup')
    BOWTIE_DIR=os.path.join(BASE_DIR, 'bowtie')
    FEATURE_DIR=os.path.join(BASE_DIR, 'feature')
    COV_DIR=os.path.join(BASE_DIR, 'bedgraph/coverage')
    NORM_COV_DIR=os.path.join(BASE_DIR, 'bedgraph/norm_coverage')
    BAM_DIR=os.path.join(BASE_DIR, 'bam')


    os.makedirs(DEDUP_DIR, exist_ok=True)
    os.makedirs(BOWTIE_DIR, exist_ok=True)
    os.makedirs(FEATURE_DIR, exist_ok=True)
    os.makedirs(COV_DIR, exist_ok=True)
    os.makedirs(NORM_COV_DIR, exist_ok=True)
    os.makedirs(BAM_DIR, exist_ok=True)
    
    futures_result=[]
    for sample in SAMPLES:

        f7=bowtie.remote(
            os.path.join(DEDUP_DIR, f"{sample}.fastq"),
            os.path.join(BOWTIE_DIR, f"{sample}.txt")
        )


        f9=feature.remote(
            os.path.join(BOWTIE_DIR, f"{sample}.txt"),
            os.path.join(FEATURE_DIR, f"{sample}.txt"),
            f7
        )
        futures_result.append(f9)

        f10=sam2bedgraph.remote(
            os.path.join(BOWTIE_DIR, f"{sample}.txt"),
            COV_DIR,
            NORM_COV_DIR,
            f7,
            ribo=True,
            offset=14
        )
        futures_result.append(f10)

        f11=sam2bam.remote(
            os.path.join(BOWTIE_DIR, f"{sample}.txt"),
            os.path.join(BAM_DIR, f"{sample}.bam"),
            f7
        )
        futures_result.append(f11)
    
    ray.get(futures_result)
    ray.shutdown()