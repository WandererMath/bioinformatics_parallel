import pandas as pd
import os

from matplotlib_venn import venn3
import matplotlib.pyplot as plt
from mylib.gtf import GTF

GTF_PATH='../reference/genomic.gtf'
gtf= GTF(GTF_PATH)


# E: phage BRU
# F: phage OCHO

# path1='/fs/ess/PAS2967/dengyw144/S21_Five/RNA-seq/deseq2/E.csv'
# path2='/fs/ess/PAS2967/dengyw144/S21_Five/RNA-seq/deseq2/F.csv'

pathC='/fs/ess/PAS2967/dengyw144/S21_Five/RNA-seq/deseq2_protein_coding/C.csv'
pathE='/fs/ess/PAS2967/dengyw144/S21_Five/RNA-seq/deseq2_protein_coding/E.csv'
pathF='/fs/ess/PAS2967/dengyw144/S21_Five/RNA-seq/deseq2_protein_coding/F.csv'
DESEQ2_PATH='../RNA-seq/deseq2_protein_coding'
df1 = pd.read_csv(pathC, sep=",", index_col=0)
df2=pd.read_csv(pathE, sep=",", index_col=0)
df3=pd.read_csv(pathF, sep=",", index_col=0)

def get_significant_genes(df, t_lfc=1, padj=0.05):
    up = df[(df['padj'] <= padj) & (df['log2FoldChange'] >= t_lfc)].index.tolist()
    down = df[(df['padj'] <= padj) & (df['log2FoldChange'] <= -t_lfc)].index.tolist()
    return up, down



up1, down1 = get_significant_genes(df1)
up2, down2 = get_significant_genes(df2)
up3, down3 = get_significant_genes(df3)


with open(os.path.join(DESEQ2_PATH, 'AND_CEF2f.txt'), 'w') as f:
    f.write('up:\n')
    
    def write_id_product(gene_list):
        tmp=[f"{g}\t{gtf.id2protein(g)}" for g in gene_list]
        f.write('\n'.join(tmp))
    write_id_product(set(up2) & set(up3) & set(up1))
    f.write('\ndown:\n')
    write_id_product(set(down2) & set(down3) & set(down1))
    f.write('\n\n')

quit()

venn3([set(up1), set(up2), set(up3)], ('C', 'E', 'F'))
plt.savefig(os.path.join(DESEQ2_PATH, 'CEF_up.pdf'))
plt.clf()

venn3([set(down1), set(down2), set(down3)], ('C', 'E', 'F'))
plt.savefig(os.path.join(DESEQ2_PATH, 'CEF_down.pdf'))
plt.clf()
quit()

with open(os.path.join(DESEQ2_PATH, 'CEF2f.txt'), 'w') as f:
    f.write('C:\n')
    
    def write_id_product(gene_list):
        tmp=[f"{g}\t{gtf.id2protein(g)}" for g in gene_list]
        f.write('\n'.join(tmp))

    f.write(f'Upregulated: {len(up1)}\n')
    write_id_product(up1)
    f.write(f'\nDownregulated: {len(down1)}\n')
    write_id_product(down1)
    f.write('\n\n')

    f.write('E:\n')
    f.write(f'Upregulated: {len(up2)}\n')
    write_id_product(up2)
    f.write(f'\nDownregulated: {len(down2)}\n')
    write_id_product(down2)
    f.write('\n\n')

    f.write('F:\n')
    f.write(f'Upregulated: {len(up3)}\n')
    write_id_product(up3)
    f.write(f'\nDownregulated: {len(down3)}\n')
    write_id_product(down3)
    f.write('\n\n')
