from mylib.bedgraph import sam2bedgraph, bedgraph_for_all_samples
from mylib.metagene import metagene_analysis_all
from mylib.gtf import GTF
PATH_SAM='../Ribo-seq/intermediate/bowtie_filtered'

# bedgraph_for_all_samples(PATH_SAM, ribo=True)
metagene_analysis_all(GTF('../reference/genomic.gtf'), '../Ribo-seq/intermediate/coverage', '../Ribo-seq/metagene')
