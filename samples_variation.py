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
from skbio.stats.distance import DistanceMatrix, permanova
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import pdist, squareform
from statsmodels.multivariate.manova import MANOVA
from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu, kruskal
 

wdir='PATH_TO_YOUR_WOKING_DIRECTORY'

Data=pd.read_excel('/'.join([wdir,'results/combined_AllData_FITvsStool_bactonly.xlsx']))
Samples=pd.read_excel('/'.join([wdir,'results/combined_protein_summary_FITvsStool.xlsx']))
annot=pd.read_csv('/'.join([wdir.removesuffix('BSA_removal'),'database/eggnog_annot.emapper.annotations']),sep='\t',skiprows=4)
annot=annot.set_index('#query')

# =============================================================================
# Format datasets
# =============================================================================

meta=Samples.loc[Samples['BSA_removal'].isin(['Stool_noBSA', 'Spin'])]
meta=meta[['SampleID','Type','Pid']].drop_duplicates().set_index('SampleID')
meta.to_csv('/'.join([wdir,'results/meta_FITvsStool.csv']), index=True)

cols=meta.index.tolist()
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

#Run Manova analysis for PCA
#add info on replicates and read in the table
PCs=pd.read_csv('/'.join([wdir,'results/PCA_FITvsStool_replicatesadded.csv']),sep=';')

group_cols = ["Type", "Pid", "Replicate"]  

for g in group_cols:
    manova_df = PCs[[g, "PC1", "PC2"]].dropna()
    maov = MANOVA.from_formula(f"PC1 + PC2 ~ {g}", data=manova_df)
    print(f'Testing for {g}:\n')
    print(maov.mv_test())
    print('\n')

# =============================================================================
# Find which PCs give individually significant differences between sample types
# and what are their loadings
# =============================================================================

def test_by_group(data, group_col):
    
    pc_cols = [c for c in data.columns if c.startswith("PC")]
    groups = data[group_col].unique()
    
    results = []

    if len(groups)==2:
        g1 = data[data[group_col] == groups[0]]
        g2 = data[data[group_col] == groups[1]]
 
        for pc in pc_cols:
            stat, pval = mannwhitneyu(g1[pc].dropna(), g2[pc].dropna(), alternative="two-sided")
            results.append({"PC": pc, "test": 'Mann-Whitney', "Stat": stat, "p_value": pval, "Group": group_col})
    
    elif len(groups)>2:
        for pc in pc_cols:
            samples = [data[data[group_col] == g][pc].dropna()for g in groups]
            
            # keep only non-empty groups
            samples = [s for s in samples if len(s) > 0]
            
            if len(samples) > 1:
                stat, pval = kruskal(*samples)
            else:
                stat, pval = None, None
            
            results.append({"PC": pc,"test": "Kruskal-Wallis","Stat": stat,"p_value": pval,"Group": group_col})
    
    results_df = pd.DataFrame(results)
    
    # FDR correction across all PCs with valid p-values
    valid = results_df["p_value"].notna()
    if valid.any():
        results_df.loc[valid, "p_fdr"] = multipletests(results_df.loc[valid, "p_value"],method="fdr_bh")[1]
    else:
        results_df["p_fdr"] = pd.NA
    
    return results_df

SignPCA = pd.DataFrame()

for group in group_cols:
    
    # Mann-Whitney test
    mw = test_by_group(PCs, group_col=group)
    SignPCA = pd.concat([SignPCA, mw], ignore_index=True)

variance=pd.DataFrame(vr,index=[f'PC{v}' for v in range(1,len(vr)+1,1)], columns=['Variance']).reset_index().rename(columns={'index':'PC'})
SignPCA=SignPCA.merge(variance, on='PC')

SignPCA.to_csv('/'.join([wdir,'results/Significance_PCA_PCs.csv']), index=False)

#Check their loadings     
loadings = pd.DataFrame(adata.varm["PCs"], index=adata.var_names, columns=[f"PC{i+1}" for i in range(variance.shape[0])])

def plot_interest(pc_int,numprot,loadings=loadings, annot=annot):
    subset=loadings[[pc_int]]
    subset['abs']=subset[pc_int].abs()
    subset=subset.merge(annot[['COG_category','Description']],left_index=True,right_index=True, how='left').reset_index()
    subset['COG_category']=subset['COG_category'].str.replace('-','unclass')
    subset['COG_category']=subset['COG_category'].fillna('unclass')

    subset['label']=subset.apply(lambda row: f'{row.COG_category};{row.Protein}', axis=1)
    subset=subset.set_index('label')
    subset=subset.sort_values(by='abs',ascending=False).head(numprot).sort_values(by=['COG_category',pc_int],ascending=False)

    
    s=subset[pc_int]
    colors = ["#D85A30" if v < 0 else "#378ADD" for v in s]
    fig, ax = plt.subplots(figsize=(numprot * 0.35 + 1,6))
    ax.vlines(s.index, 0, s.values, colors=colors, linewidth=2)
    ax.scatter(s.index, s.values,  color=colors, zorder=3, s=30)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel(f"{pc_int} loading")
    ax.set_title(f"Top {numprot} protein loadings in {pc_int}")
    ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
    sb.despine()
    plt.tight_layout()
    plt.show()
    
    return subset[[pc_int,'COG_category']]

pc3top=plot_interest('PC3',100)
pc3top.to_csv('/'.join([wdir,'results/PCA_FITvsStool_PC3_loadings_top100.csv']), index=True)

pc3top['COG_category'].value_counts()

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
# Run PERMANOVA analysis
# =============================================================================

group_cols = ["Type", "Pid"]  

#make sure indexes are in the same order
meta=meta.loc[DataNorm.index.tolist()]

Permanova=pd.DataFrame()
for g in group_cols:
    permanova_res = permanova(distance_matrix=dm, grouping=meta[g],permutations=999)
    
    F=permanova_res["test statistic"]
    between=permanova_res["number of groups"] - 1
    within=permanova_res["sample size"] - permanova_res["number of groups"]

    # Compute R²
    R2=(F*between)/(F*between+within)
    
    perm=pd.DataFrame({'Group':[g],'PseudoF':[F], 'R2':[R2], 'Pval':[permanova_res["p-value"]]})
    Permanova=pd.concat([Permanova, perm])
    
Permanova.to_csv('/'.join([wdir,'results/Permanova_Jaccard_FITvsStool.csv']), index=False)