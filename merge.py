import os

from mylib.utils import get_paths_ends_with_something
from func import merge_bedgraphs

for dir_path, prefix in zip(['RNA-seq/bam', 'Ribo-seq/bam'], ["RNA", "Ribo"]):
    args={v:[v+j for j in '123'] for v in 'ABC'}

    for sample, replicates in args.items():
        for strand in ['plus', 'minus']:
            paths=[os.path.join(dir_path, f'{replicate}.{strand}.bedgraph') for replicate in replicates]
            out_path=os.path.join(dir_path, f'{prefix}-{sample}.{strand}.bedgraph')
            merge_bedgraphs(paths, out_path)


