for f in *.bam; do
    out_plus="${f%.bam}.plus.bedgraph"
    bedtools genomecov -ibam $f -bg -strand + > $out_plus
    out_minus="${f%.bam}.minus.bedgraph"
    bedtools genomecov -ibam $f -bg -strand - > $out_minus
done