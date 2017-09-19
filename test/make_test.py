import pandas as pd

#df = pd.read_table('/stor/work/Lambowitz/cdw2854/miRNA/tgirt_data/bam_files/miRNA_counts.tsv') \
df = pd.read_table('/stor/work/Lambowitz/cdw2854/miRNA/tgirt_data/hand_mix_new_count.tsv') \
        .assign(cpm = lambda d: (d.HM_NTT1 + d.HM_NTT2)/2)\
        .rename(columns={'Row.names':'ID'})\
        .pipe(lambda d: d[['ID','cpm']])
df.to_csv('tgirt_count.csv',header=False, index=False)
