import os
import multiprocessing
import subprocess
import fcntl
import ray
import pysam
from mylib.utils import get_paths_ends_with_something


N_CPU=multiprocessing.cpu_count()

BOWTIE_REF_PATH='reference/ecoli'
GTF_PATH='reference/genomic.gtf'
MIN_INS=0
MAX_INS=100
@ray.remote
def unzip(path, output_path):
    command=f"gunzip -c {path} > {output_path}"
    subprocess.run([command], shell=True, check=True)
    return 0

@ray.remote(num_cpus=12)
def trim(path_r1, path_r2, output_path_r1, output_path_r2, *_):
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
def merge(path_r1, path_r2, output_name, *_):
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
def move_merged(path_name, output_folder, *_):
    command=f'mv {path_name}.assembled.fastq {output_folder}'
    subprocess.run([command], shell=True, check=True)
    return 0

@ray.remote(num_cpus=12)
def dedup(path_name, output_path, *_):
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
def bowtie(path_name, output_path, *_):
    command = "bowtie" + \
            " -m 1"+\
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
def filter_1_sam(path, output_path, low_threshold, high_threshold, *_):
    # Open SAM input and output
    with pysam.AlignmentFile(path, "r") as infile, \
        pysam.AlignmentFile(output_path, "w", template=infile) as outfile:
        for read in infile:
            if not read.is_unmapped and read.query_length >= low_threshold\
                and read.query_length<=high_threshold:
                outfile.write(read)         

@ray.remote(num_cpus=12)
def feature(path, output_path, *_):
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

def write_log(command, r):
    with open('out.log', 'a') as f: 
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(command+'\n')
        if r.stdout is not None:
            f.write(r.stdout)
        if r.stderr is not None:
            f.write(r.stderr)
        f.write('\n')
        fcntl.flock(f, fcntl.LOCK_UN)

@ray.remote
def sam2bam(sam_path, bam_path, *_):
    cmd=f'samtools view -S -b {sam_path} | samtools sort -o {bam_path}'
    cmd2=f'samtools index {bam_path}'
    r=subprocess.run([cmd], shell=True, check=True,  capture_output=True, text=True)
    write_log(cmd, r)
    r2=subprocess.run([cmd2], shell=True, check=True,  capture_output=True, text=True)
    write_log(cmd2, r2)


@ray.remote
def sam2bedgraph(sam_path, out_dir, norm_out_dir, *_, ribo=True, offset=14):
    sam_path=os.path.abspath(sam_path)
    FIRST_LINE='track type=bedGraph\n'
    BASENAME=os.path.basename(sam_path).split('.')[0]

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(norm_out_dir, exist_ok=True)
    prefix='Ribo-' if ribo else 'RNA-'

    OUTPUT_PLUS=os.path.join(out_dir, prefix+BASENAME + '-plus.bedgraph')
    OUTPUT_MINUS=os.path.join(out_dir, prefix+BASENAME + '-minus.bedgraph')

    OUTPUT_PLUS_NORM=os.path.join(norm_out_dir, prefix+BASENAME + '-plus.bedgraph')
    OUTPUT_MINUS_NORM=os.path.join(norm_out_dir, prefix+BASENAME + '-minus.bedgraph')

    result={'+':{}, '-':{}} # chrom: dict
                #     # dict :=   start: value
    with pysam.AlignmentFile(sam_path, "r") as samfile:
        for read in samfile.fetch():
            chrom, start, end = read.reference_name, read.reference_start, read.reference_end
            
            strand = '-' if read.is_reverse else '+'   
            if read.is_unmapped:
                continue
            if chrom not in result[strand]:
                result[strand][chrom] = {} 
                
            if ribo:
                start+=1
                if strand=='+':
                    pos= end-offset
                else:
                    pos= start+offset 
            else:
                pos= (start + end)//2

            if pos not in result[strand][chrom]:
                result[strand][chrom][pos] = 1
            else:
                result[strand][chrom][pos] += 1

    total=0
    with open(OUTPUT_PLUS, 'w') as f:
        f.write(FIRST_LINE)
        for chrom in result['+']:
            for pos in result['+'][chrom]:
                total+=result['+'][chrom][pos]
                f.write(f"{chrom}\t{pos}\t{pos+1}\t{result['+'][chrom][pos]}\n")

    with open(OUTPUT_MINUS, 'w') as f:
        f.write(FIRST_LINE)
        for chrom in result['-']:
            for pos in result['-'][chrom]:
                total+=result['-'][chrom][pos]
                f.write(f"{chrom}\t{pos}\t{pos+1}\t{result['-'][chrom][pos]}\n")
    
    # Normalize 
    total/=(1E6)
    with open(OUTPUT_PLUS_NORM, 'w') as f:
        f.write(FIRST_LINE)
        for chrom in result['+']:
            for pos in result['+'][chrom]:
                result['+'][chrom][pos]/=total
                f.write(f"{chrom}\t{pos}\t{pos+1}\t{result['+'][chrom][pos]}\n")

    with open(OUTPUT_MINUS_NORM, 'w') as f:
        f.write(FIRST_LINE)
        for chrom in result['-']:
            for pos in result['-'][chrom]:
                result['-'][chrom][pos]/=total
                f.write(f"{chrom}\t{pos}\t{pos+1}\t{result['-'][chrom][pos]}\n")
    # breakpoint()
    # return result


def merge_bedgraphs_custom(input_samples_list, output_path):
    result={}
    for path in input_samples_list:
        with open(path) as f:
            next(f)
            for line in f:
                try:
                    x, y, z, value= line.strip().split('\t')
                    value=float(value)
                    if (x, y, z) not in result:
                        result[(x, y, z)]=0
                    result[(x, y, z)]+=value
                except:
                    continue
    with open(output_path, 'w') as f:
        f.write('track type=bedGraph\n')
        for key in result:
            f.write(f"{key[0]}\t{key[1]}\t{key[2]}\t{result[key]}\n")


def merge_bedgraphs(input_samples_list, output_path):
    cmd=f'bedtools unionbedg -i {' '.join(input_samples_list)} > {output_path}.tmp'
    cmd_sum="awk '{print $1, $2, $3, $4+$5+$6}' OFS='\t' "+f"{output_path}.tmp > {output_path}"
    r=subprocess.run([cmd], shell=True, check=True,  capture_output=True, text=True)
    r2=subprocess.run([cmd_sum], shell=True, check=True,  capture_output=True, text=True)
    write_log(cmd, r)
    write_log(cmd_sum, r2)
    