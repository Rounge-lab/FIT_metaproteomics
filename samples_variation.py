#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 11:43:33 2026

PCA analysis and Jaccard distance/PCoA between samples

@author: ekateria
"""
import pandas as pd
import seaborn as sb
import scanpy as sc
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import pdist, squareform
from itertools import combinations as comb
from statannotations.Annotator import Annotator as ann 

wdir='PATH_TO_YOUR_WOKING_DIRECTORY'

Data=pd.read_excel(f'{wdir}/input/true_MaxLFQ_nohuman.xlsx')
Samples=pd.read_excel(f'{wdir}/results/combined_protein_summary.xlsx')
meta=pd.read_excel(f'{wdir}/input/meta.xlsx')
annot=pd.read_csv(f'{wdir}/database/eggnog_annot.emapper.annotations',sep='\t',skiprows=4)
annot=annot.set_index('#query')

# =============================================================================
# Format datasets
# =============================================================================

cols=meta['SampleID']
cols.append('Protein')

DataNorm=Data[cols].set_index('Protein').T

#keep proteins that are detected at least in one sample 
Prev=DataNorm.copy()
Prev[Prev>0]=1
Sums=Prev.sum(axis=0)
tokeep=Sums[Sums>0].index.tolist()
DataNorm=DataNorm[tokeep]

#arrange the samples in the same order as meta
DataNorm=DataNorm.loc[meta.index]

#Create annotated data
adata = ad.AnnData(X=DataNorm.values.astype(float), obs=meta.copy(), var=pd.DataFrame(index=DataNorm.columns))

# =============================================================================
# Impute zeros for PCA: 0.5 of minimum intensity value for this protein; 
# Log-transform & Z-scale the data 
# =============================================================================

X=adata.X.copy()
X[X<=0]=np.nan
protein_mins=np.nanmin(X,axis=0)
X_imp=np.where(np.isnan(X),0.5*protein_mins,X)
adata.X=X_imp

#Log2-transform the data
adata.X = np.log2(adata.X)

#Median-center each sample to account for multiple runs
sam_medians=np.median(adata.X,axis=1,keepdims=True)
adata.X=adata.X-sam_medians

#Z-score scale the proteins
sc.pp.scale(adata,zero_center=True,max_value=10)

# =============================================================================
# Run PCA and make a plot
# =============================================================================

sc.tl.pca(adata)
#Get explained variance
vr = adata.uns['pca']['variance_ratio'] * 100

PCvals = adata.obsm['X_pca'][:, :2]
PCs=adata.obs[['Pid', 'Type']].copy()

PCvals=adata.obsm['X_pca'][:,:len(vr)+1]
for v in range(1,len(vr)+1,1):
    PCs[f'PC{v}']=PCvals[:,v-1]

PCs.to_csv('/'.join([wdir,'results/PCA_FITvsStool.csv']), index=True)

fig=sb.scatterplot(PCs,x='PC1',y='PC2',hue='Type',style='Pid', palette=['#7ebfaa','#ec9e7d'],s=180)

xmin, xmax = fig.get_xlim()
ymin, ymax = fig.get_ylim()

fig.vlines(x=0, ymin=ymin, ymax=ymax, color='grey', linestyle='--', linewidth=1)
fig.hlines(y=0, xmin=xmin, xmax=xmax, colors='grey', linestyle='--', linewidth=1)

fig.set(xlabel=f'PC1 ({vr[0]:.1f}%)', ylabel=f'PC2 ({vr[1]:.1f}%)')

plt.savefig('/'.join([wdir,'results/PCA_FITvsStool.pdf']))

# =============================================================================
# Calculate Jaccard index and run PCoA
# =============================================================================
Jac=squareform(pdist(Prev.values,metric='jaccard'))
dm=DistanceMatrix(Jac,ids=DataNorm.index.astype(str))

jac_ordin=pcoa(dm)

pcoa_df = jac_ordin.samples.iloc[:, :2].copy()
pcoa_df.columns = ["PCo1", "PCo2"]

# Variance explained
ve = (jac_ordin.proportion_explained.iloc[:2] * 100).values

# Merge with metadata for plotting
plot_df = meta.loc[pcoa_df.index].copy()
plot_df = plot_df.join(pcoa_df)

fig=sb.scatterplot(plot_df,x='PCo1',y='PCo2',hue='Type',style='Pid', palette=['#7ebfaa','#ec9e7d'],s=180)

xmin, xmax = fig.get_xlim()
ymin, ymax = fig.get_ylim()

fig.vlines(x=0, ymin=ymin, ymax=ymax, color='grey', linestyle='--', linewidth=1)
fig.hlines(y=0, xmin=xmin, xmax=xmax, colors='grey', linestyle='--', linewidth=1)

fig.set(xlabel=f'PCo1 ({ve[0]:.1f}%)', ylabel=f'PCo2 ({ve[1]:.1f}%)')
plt.savefig('/'.join([wdir,'results/PCoA_Jaccard_ITvsStool.pdf']))

Jaccard=pd.DataFrame(Jac, columns=Prev.index.tolist(), index=Prev.index.tolist())
Jaccard.to_csv('/'.join([wdir,'results/Jaccard_FITvsStool.csv']), index=True)

# =============================================================================
# Find an average Jacacrd distance between ID to ID; types and replicates
# =============================================================================

Jac_pairs=pd.DataFrame()
meta=meta.reset_index()
#Find Jaccard pairs within each ID
allcomp=list(comb(meta['SampleID'].tolist(),2))
for a in allcomp:
    if meta.loc[meta['SampleID']==a[0],'Pid'].item()==meta.loc[meta['SampleID']==a[1],'Pid'].item():
        idcat='replicates'
    else:
        idcat='ID to ID'
    jac={'S1':[a[0]],'S2':[a[1]],'jac':[Jaccard.loc[a[0],a[1]]], 
         'Type1':[meta.loc[meta['SampleID']==a[0],'Type'].item()],'Type2':[meta.loc[meta['SampleID']==a[1],'Type'].item()],
         'IDcat':[idcat]}
    Jac_pairs=pd.concat([Jac_pairs, pd.DataFrame(jac)])
    
Jac_pairs['Typecat']=Jac_pairs['Type1']+'_vs_'+Jac_pairs['Type2']
Jac_pairs['Typecat']=Jac_pairs['Typecat'].apply(lambda row: 'FIT_vs_Stool' if row=='Stool_vs_FIT' else row)

palette={'replicates':'#9697b1','ID to ID':'#acc0c3'}
x_order=['FIT_vs_FIT','Stool_vs_Stool','FIT_vs_Stool']
hue_order=['replicates','ID to ID']

pairs = []

for p in x_order:
    for c1, c2 in comb(hue_order, 2):
        pairs.append(((p, c1), (p, c2)))

print(pairs)

ax=sb.boxplot(Jac_pairs, y='Typecat',x='jac', hue='IDcat',palette=palette,hue_order=hue_order,order=x_order,legend=False)

annj=ann(ax,pairs,data=Jac_pairs,y='Typecat',x='jac', hue='IDcat', order=x_order, hue_order=hue_order,orient='h')
annj.configure(test='Mann-Whitney',text_format='star',loc='inside',hide_non_significant=True)
annj.apply_and_annotate()   
ax.set(xlabel='Jaccard distance',ylabel='')
ax.set_yticklabels(['FIT','Stool','FIT vs Stool'])
plt.savefig(f'{wdir}/results/Jaccard_by_Type_and_replicates.pdf')