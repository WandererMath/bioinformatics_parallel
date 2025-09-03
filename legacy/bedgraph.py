from mylib.bedgraph import sam2bedgraph, bedgraph_for_all_samples


PATH_SAM='../RNA-seq/intermediate/bowtie_filtered'

bedgraph_for_all_samples(PATH_SAM, ribo=False)

