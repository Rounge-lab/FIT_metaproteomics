# -*- coding: utf-8 -*-
"""

Check taxonomy of the detected genes

"""

import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator as ann


wdir='PATH_TO_YOUR_WORKING_DIRECTORY'

Final=pd.read_excel(f'{wdir}/input/true_MaxLFQ_nohuman.xlsx')
meta=pd.read_csv(f'{wdir}/input/meta.csv',sep=',')
Samples=pd.read_excel(f'{wdir}/results/combined_protein_summary.xlsx')


Final=Final.set_index('Protein')

tax_file=pd.read_csv(f'{wdir}/Gene_MAG_list.csv', sep=',')
tax_file=tax_file.set_index('GeneIDs')

#Keep only FIT-spinned and stool samples
final_cols=Final.columns.tolist()
final_cols=[f for f in final_cols if f in meta['SampleID'].tolist()]
Final=Final[final_cols]

#Make taxonomy presence/absence data
Taxa=Final.copy()
Taxa[Taxa>0]=1
Taxa=Taxa.merge(tax_file, left_index=True, right_index=True, how='left')

Taxa['phylum']=Taxa['classification'].apply(lambda row: row.split('p__')[1])
Taxa['phylum']=Taxa['phylum'].apply(lambda row: row.split(';')[0])
Taxa['phylum']=Taxa['phylum'].apply(lambda row: row.split('_')[0])

#Rename phyla to comply with updated nomenclature
Taxa.loc[Taxa['phylum']=='Firmicutes','phylum']='Bacillota'
Taxa.loc[Taxa['phylum']=='Proteobacteria','phylum']='Pseudomonadota'
Taxa.loc[Taxa['phylum']=='Actinobacteriota','phylum']='Actinomycetota'
Taxa.loc[Taxa['phylum']=='ABacteroidetes','phylum']='Bacteroidota'

Taxa['species']=Taxa['classification'].apply(lambda row: row.split('s__')[-1])
len(Taxa['species'].dropna().unique().tolist())

#Make a summary of which phyla are detected in each species, and how many species belonging to that phylum

#Make per-phylum-per-species summary
Taxa['species']=Taxa['species'].apply(lambda row: np.nan if len(row)<2 else row)
Taxa_knownsp=Taxa.dropna(subset='species')

Taxa_melt=Taxa_knownsp.melt(id_vars=['phylum','species'],value_vars=final_cols,value_name='presence',var_name='SampleID')
Taxa_melt=Taxa_melt.query('presence==1').drop(columns='presence')

Ph_pres=Taxa_melt.groupby(['SampleID','phylum'])['species'].nunique().reset_index(name='NumSpecies')
Ph_pres=Ph_pres.merge(meta,on='SampleID', how='left')

#Calculate phylum prevalence for each sample type
Ph_preval=Ph_pres.groupby('Type')['phylum'].value_counts().reset_index()
Ph_preval.loc[Ph_preval['Type']=='FIT','Preval']=Ph_preval.loc[Ph_preval['Type']=='FIT','count']/len(meta.loc[meta['Type']=='FIT'])*100
Ph_preval.loc[Ph_preval['Type']=='Stool','Preval']=Ph_preval.loc[Ph_preval['Type']=='Stool','count']/len(meta.loc[meta['Type']=='Stool'])*100

order=['Actinomycetota','Bacteroidota','Bacillota','Pseudomonadota','Desulfobacterota','Verrucomicrobiota','Fusobacteriota','Cyanobacteria','Synergistota','Campylobacterota','Thermoplasmatota']

palette={'Stool':'#7ebfaa','FIT':'#ec9e7d'}
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(8,6),sharey=True)
sb.barplot(Ph_preval,x='Preval',y='phylum',hue='Type',hue_order=['FIT','Stool'],palette=palette,ax=ax1,order=order,legend=False)
ax1.set(ylabel='',xlabel='Prevalence, %')

sb.boxplot(Ph_pres,x='NumSpecies',y='phylum',hue='Type',hue_order=['FIT','Stool'],palette=palette,log_scale=True,legend=False,ax=ax2,order=order)
ax2.grid(True, axis='x', linestyle='--', alpha=0.5)
ax2.set(ylabel='',xlabel='Number of species')
plt.tight_layout()

pairs = []

for p in order:
    pairs.append(((p, 'FIT'), (p, 'Stool')))

print(pairs)

annj=ann(ax2,pairs,data=Ph_pres,x='NumSpecies',y='phylum',hue='Type',hue_order=['FIT','Stool'],order=order,orient='h')
annj.configure(test='Mann-Whitney',text_format='star',loc='inside',hide_non_significant=True)
annj.apply_and_annotate()  

plt.savefig(f'{wdir}/results/Taxonomy_annot.pdf')
