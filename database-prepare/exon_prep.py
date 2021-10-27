#%%
from autocnv import settings
from autocnv.database import DataBase
from gtfparse import read_gtf
import pandas as pd
refGene_file = '/Users/zhonghua/data/gene-toolkit/raw-data/hg19.refGene.gtf.gz'
#%%
refGene_df = read_gtf(refGene_file)
refGene_df = refGene_df[refGene_df['feature'] == 'exon']
#%%
gene_db = pd.read_csv(settings.GENE_DATABASE, sep='\t')
#%%

merge_df = gene_db.merge(refGene_df, left_on='transcript', right_on='transcript_id')[[
    '#chrom', 'start_y', 'end_y', 'gene_id_x', 'symbol', 'exon_number', 'transcript'
]]
merge_df.columns = [x.strip('_x').strip('_y') for x in merge_df.columns]

merge_df.to_csv(settings.GENE_EXON_DATABASE.strip('.gz'), sep='\t', index=False)
#%%


