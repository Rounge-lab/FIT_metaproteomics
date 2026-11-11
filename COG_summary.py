# -*- coding: utf-8 -*-
"""

Make COG plot by samples

"""
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator as annot

wdir='PATH_TO_YOUR_WORKING_DIRECTORY'

dbannot=pd.read_csv(f'{wdir}/database/eggnog_annot.emapper.annotations.csv',sep='\t',skiprows=4)
dbannot=dbannot.set_index('#query')

Final=pd.read_excel(f'{wdir}/input/true_MaxLFQ_nohuman.xlsx')
meta=pd.read_excel(f'{wdir}/input/meta.xlsx')

Final=Final.set_index('Protein')
Final[Final>0]=1
COGs=Final.merge(dbannot[['COG_category']],left_index=True,right_index=True, how='left')

COG_melt=COGs.melt(id_vars=['COG_category'],value_vars=final_cols,value_name='presence',var_name='SampleID')
COG_melt=COG_melt.query('presence==1').drop(columns='presence')

COG_pres=COG_melt.groupby(['SampleID','COG_category']).value_counts().reset_index(name='Num_Proteins')
COG_pres['Proportion'] = (COG_pres['Num_Proteins'] / COG_pres.groupby('SampleID')['Num_Proteins'].transform('sum')* 100)
COG_pres=COG_pres.merge(meta,on='SampleID', how='left')

#Find total number of proteins in each category
numprot=COGs['COG_category'].value_counts().reset_index()

#Plot only top10 prevalent COGs
COG_pres=COG_pres.merge(numprot,on='COG_category',how='left')
COG_pres=COG_pres.rename(columns={'count':'Total'})
COG_pres['COG']=COG_pres.apply(lambda row: f'{row.COG_category} (Nprot={row.Total})',axis=1)
order = COG_pres.groupby('COG')['Proportion'].median().sort_values(ascending=False).reset_index()


palette={'Stool':'#7ebfaa','FIT':'#ec9e7d'}

plt.subplots(figsize=(8,4))
ax=sb.boxplot(COG_pres.loc[COG_pres['COG'].isin(order.loc[0:10,'COG'].tolist())],x='COG',y='Proportion',hue='Type',palette=palette, legend=False,hue_order=['FIT','Stool'],order=order.loc[0:10,'COG'].tolist())
ax.set_xticklabels(ax.get_xticklabels(),fontsize=8,rotation=90)
ax.set_yticklabels(ax.get_yticklabels(),fontsize=8)
ax.set(xlabel='')
ax.set_ylabel('Proportion,%',fontsize=8)
#Make separation lines between categories
n = len(ax.get_xticks())
for x in range(n - 1):
    ax.axvline(x + 0.5, linestyle='--', color='gray', linewidth=0.8, alpha=0.7)

pairs = []

for c in order.loc[0:10,'COG'].tolist():
    pairs.append(((c, 'FIT'), (c, 'Stool')))

ann=annot(ax,pairs,data=COG_pres.loc[COG_pres['COG'].isin(order.loc[0:10,'COG'].tolist())],x='COG',y='Proportion', hue='Type', hue_oder=['FIT','Stool'])
ann.configure(test='Mann-Whitney',text_format='star',loc='inside',hide_non_significant=True)

ax, annotations = ann.apply_and_annotate()

plt.savefig(f'{wdir}/results/Fraction_by_COGs_annot.pdf')

# extract stats into tabular form
rows = []
for ann in annotations:
    res = ann.data  # StatResult

    rows.append({
        "group1": ann.structs[0]["label"],
        "group2": ann.structs[1]["label"],
        "test": getattr(res, "test_description", None),
        "pvalue": getattr(res, "pvalue", None),
    })

MannWStats = pd.DataFrame(rows)

MannWStats.to_csv(f'{wdir}/results/Fraction_by_COG_MannW.csv',sep='\t', index=False)
