import os
import multiprocessing
import subprocess
import fcntl
import ray
from mylib.utils import get_paths_ends_with_something

N_CPU=multiprocessing.cpu_count()
BOWTIE_REF_PATH='reference/MyStrain'
MIN_INS=23
MAX_INS=32
@ray.remote
def unzip(path, output_path):
    command=f"gunzip -c {path} > {output_path}"
    subprocess.run([command], shell=True, check=True)
    return 0

@ray.remote(num_cpus=12)
def trim(path_r1, path_r2, output_path_r1, output_path_r2, *ready):
    command=f'cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o {output_path_r1} \
        -p {output_path_r2} \
        --discard-untrimmed --cores=12 --pair-filter=any --minimum-length=20 --maximum-length=55 -q 20,20 \
        {path_r1}\
        {path_r2}'
    r=subprocess.run([command], shell=True,  capture_output=True, text=True)
    with open('out.log', 'a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(command+'\n')
        f.write(r.stdout)
        f.write(r.stderr)
        f.write('\n')
        fcntl.flock(f, fcntl.LOCK_UN)
        
        
    return 0


@ray.remote(num_cpus=12)
def merge(path_r1, path_r2, output_name, *ready):
    command=f'pear -f {path_r1} \
        -r {path_r2}\
        -o {output_name} \
        -n 20 --threads 12 --memory 4G'
    r=subprocess.run([command], shell=True, check=True,  capture_output=True, text=True)

    with open('out.log', 'a') as f:    
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(command+'\n')
        f.write(r.stdout)
        f.write(r.stderr)
        f.write('\n')
        fcntl.flock(f, fcntl.LOCK_UN)
    return 0

@ray.remote
def move_merged(path_name, output_folder, *ready):
    command=f'mv {path_name}.assembled.fastq {output_folder}'
    subprocess.run([command], shell=True, check=True)
    return 0

@ray.remote(num_cpus=12)
def dedup(path_name, output_path, *ready):
    command = "fastp" + \
            " -i " + path_name + \
            " -o " + output_path + \
            " --dedup" + \
            " -Q -A -G" + \
            " -j /dev/null -h /dev/null" + \
            " --thread 12" + \
            " --trim_poly_x" + \
            " --length_required 20"
    r= subprocess.run([command], shell=True, check=True)
    with open('out.log', 'a') as f:    
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(command+'\n')
        f.write(r.stdout)
        f.write(r.stderr)
        f.write('\n')
        fcntl.flock(f, fcntl.LOCK_UN)
    return 0

@ray.remote(num_cpus=12)
def bowtie(path_name, output_path, *ready):
    command = "bowtie" + \
            " -q" + \
            " -3 5" + \
            " -x " + BOWTIE_REF_PATH + \
            " --best -S" + \
            " --threads 12" + \
            " --minins" + MIN_INS + \
            " --maxins" + MAX_INS + \
            " " + path_name+ \
            " " + output_path

    r=subprocess.run([command], shell=True, check=True,  capture_output=True, text=True)

    with open('out.log', 'a') as f: 
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(command+'\n')
        f.write(r.stdout)
        f.write(r.stderr)
        f.write('\n')
        fcntl.flock(f, fcntl.LOCK_UN)
    return 0

BASE_DIR='RNA-seq/intermediate'
BASE_DIR='test_data'
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
    os.makedirs(UNZIPPED_DIR, exist_ok=True)
    os.makedirs(TRIMMED_DIR, exist_ok=True)
    os.makedirs(MERGED_DIR, exist_ok=True)
    os.makedirs(MERGED_ASSEMBLED_DIR, exist_ok=True)
    os.makedirs(DEDUP_DIR, exist_ok=True)
    os.makedirs(BOWTIE_DIR, exist_ok=True)
    
    futures_result=[]
    for sample in SAMPLES:
        # f1=unzip.remote(
        #     os.path.join(RAW_DIR, f"{sample}_R1.fastq.gz"),
        #     os.path.join(UNZIPPED_DIR, f"{sample}_R1.fastq")
        # )
        # f2=unzip.remote(
        #     os.path.join(RAW_DIR, f"{sample}_R2.fastq.gz"),
        #     os.path.join(UNZIPPED_DIR, f"{sample}_R2.fastq")
        # )
        # f3=trim.remote(os.path.join(UNZIPPED_DIR, f"{sample}_R1.fastq"),
        #     os.path.join(UNZIPPED_DIR, f"{sample}_R2.fastq"),
        #     os.path.join(TRIMMED_DIR, f"{sample}_R1.fastq"),
        #     os.path.join(TRIMMED_DIR, f"{sample}_R2.fastq"), f1, f2
        # )
        # f4=merge.remote(
        #     os.path.join(TRIMMED_DIR, f"{sample}_R1.fastq"),
        #     os.path.join(TRIMMED_DIR, f"{sample}_R2.fastq"),
        #     os.path.join(MERGED_DIR, sample),
        #     # f3
        # )
        # f5=move_merged.remote(
        #     os.path.join(MERGED_DIR, sample),
        #     MERGED_ASSEMBLED_DIR,
        #     f4
        # )
        f6=dedup.remote(
            os.path.join(MERGED_ASSEMBLED_DIR, f"{sample}.assembled.fastq"),
            os.path.join(DEDUP_DIR, f"{sample}.fastq"),
            # f5
        )
        f7=bowtie.remote(
            os.path.join(DEDUP_DIR, f"{sample}.fastq"),
            os.path.join(BOWTIE_DIR, f"{sample}.txt"),
            f6
        )
        futures_result.append(f7)
    
    ray.get(futures_result)
    ray.shutdown()