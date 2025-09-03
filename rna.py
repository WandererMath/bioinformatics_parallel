import os
import multiprocessing
import subprocess
import fcntl
import ray
import pysam
from mylib.utils import get_paths_ends_with_something
from func import *

BASE_DIR='RNA-seq/intermediate'
# BASE_DIR='test_data'

if __name__=='__main__':
    print('N_CPU: ', N_CPU)
    ray.init(num_cpus=N_CPU)
    
    # retrieve sample names
    SAMPLES=[]
    RAW_DIR=os.path.join(BASE_DIR, 'raw')
    RAW_FILES=get_paths_ends_with_something(RAW_DIR, '.fastq.gz')
    RAW_FILES_BASE=[os.path.basename(f) for f in RAW_FILES]
    SAMPLES=list(set([f.split('_')[0] for f in RAW_FILES_BASE]))
    print('N_SAMPLES: ', len(SAMPLES))
    print('SAMPLES: ', SAMPLES)
    
    UNZIPPED_DIR=os.path.join(BASE_DIR, 'unzipped')
    TRIMMED_DIR=os.path.join(BASE_DIR, 'trimmed')
    MERGED_DIR=os.path.join(BASE_DIR, 'merged')
    MERGED_ASSEMBLED_DIR=os.path.join(BASE_DIR, 'merged_assembled')
    DEDUP_DIR=os.path.join(BASE_DIR, 'dedup')
    BOWTIE_DIR=os.path.join(BASE_DIR, 'bowtie')
    BOWTIE_FILTERED_DIR=os.path.join(BASE_DIR, 'bowtie_filtered')
    FEATURE_DIR=os.path.join(BASE_DIR, 'feature')
    os.makedirs(UNZIPPED_DIR, exist_ok=True)
    os.makedirs(TRIMMED_DIR, exist_ok=True)
    os.makedirs(MERGED_DIR, exist_ok=True)
    os.makedirs(MERGED_ASSEMBLED_DIR, exist_ok=True)
    os.makedirs(DEDUP_DIR, exist_ok=True)
    os.makedirs(BOWTIE_DIR, exist_ok=True)
    os.makedirs(BOWTIE_FILTERED_DIR, exist_ok=True)
    os.makedirs(FEATURE_DIR, exist_ok=True)
    
    futures_result=[]
    for sample in SAMPLES:
        f1=unzip.remote(
            os.path.join(RAW_DIR, f"{sample}_R1.fastq.gz"),
            os.path.join(UNZIPPED_DIR, f"{sample}_R1.fastq")
        )
        f2=unzip.remote(
            os.path.join(RAW_DIR, f"{sample}_R2.fastq.gz"),
            os.path.join(UNZIPPED_DIR, f"{sample}_R2.fastq")
        )
        f3=trim.remote(os.path.join(UNZIPPED_DIR, f"{sample}_R1.fastq"),
            os.path.join(UNZIPPED_DIR, f"{sample}_R2.fastq"),
            os.path.join(TRIMMED_DIR, f"{sample}_R1.fastq"),
            os.path.join(TRIMMED_DIR, f"{sample}_R2.fastq"), f1, f2
        )
        f4=merge.remote(
            os.path.join(TRIMMED_DIR, f"{sample}_R1.fastq"),
            os.path.join(TRIMMED_DIR, f"{sample}_R2.fastq"),
            os.path.join(MERGED_DIR, sample),
            f3
        )
        f5=move_merged.remote(
            os.path.join(MERGED_DIR, sample),
            MERGED_ASSEMBLED_DIR,
            f4
        )
        f6=dedup.remote(
            os.path.join(MERGED_ASSEMBLED_DIR, f"{sample}.assembled.fastq"),
            os.path.join(DEDUP_DIR, f"{sample}.fastq"),
            f5
        )
        f7=bowtie.remote(
            os.path.join(DEDUP_DIR, f"{sample}.fastq"),
            os.path.join(BOWTIE_DIR, f"{sample}.txt"),
            f6
        )

        f8=filter_1_sam.remote(
            os.path.join(BOWTIE_DIR, f"{sample}.txt"),
            os.path.join(BOWTIE_FILTERED_DIR, f"{sample}.txt"),
            23,
            100,
            f7
        )
        f9=feature.remote(
            os.path.join(BOWTIE_FILTERED_DIR, f"{sample}.txt"),
            os.path.join(FEATURE_DIR, f"{sample}.txt"),
            f8
        )
        futures_result.append(f9)
    
    ray.get(futures_result)
    ray.shutdown()