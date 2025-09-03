import os
import multiprocessing
import subprocess
import fcntl
import ray
import pysam
from mylib.utils import get_paths_ends_with_something


N_CPU=multiprocessing.cpu_count()

BOWTIE_REF_PATH='reference/MyStrain'
GTF_PATH='reference/genomic.gtf'
MIN_INS=0
MAX_INS=100
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
        if r.stdout is not None:
            f.write(r.stdout)
        if r.stderr is not None:
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
        if r.stdout is not None:
            f.write(r.stdout)
        if r.stderr is not None:
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
        if r.stdout is not None:
            f.write(r.stdout)
        if r.stderr is not None:
            f.write(r.stderr)
        f.write('\n')
        fcntl.flock(f, fcntl.LOCK_UN)
    return 0

@ray.remote(num_cpus=12)
def bowtie(path_name, output_path, *ready):
    command = "bowtie" + \
            " -q" + \
            " -3 5 " + \
            BOWTIE_REF_PATH + \
            " --best -S" + \
            " --threads 12" + \
            f" --minins {MIN_INS}" + \
            f" --maxins {MAX_INS}" + \
            " " + path_name+ \
            " " + output_path

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
    return 0


@ray.remote
def filter_1_sam(path, output_path, low_threshold, high_threshold, *ready):
    # Open SAM input and output
    with pysam.AlignmentFile(path, "r") as infile, \
        pysam.AlignmentFile(output_path, "w", template=infile) as outfile:
        for read in infile:
            if not read.is_unmapped and read.query_length >= low_threshold\
                and read.query_length<=high_threshold:
                outfile.write(read)         

@ray.remote(num_cpus=12)
def feature(path, output_path, *ready):
    command=f'featureCounts  -s 1 -M -Q 20 -T 12 -t gene -g gene_id -a {GTF_PATH} -o {output_path} {path}'
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
    return 0


