#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 10:53:22 2026

Check hydrophobicity of the top100 proteins

@author: ekaterinaavershina
"""
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import seaborn as sb
import matplotlib.pyplot as plt

wdir='PATH_TO_YOUR_WORKING_DIR'
seqsfile='PATH_TO_PROTEIN_SEQUENCE_FASTA'


pc3top=pd.read_csv('/'.join([wdir,'results/PCA_FITvsStool_PC3_loadings_top100.csv']))

for ix,p in pc3top.iterrows():
    print('------------------------------------\n')
    print(f"Calculating GRAVY for {p.label.split(';')[1]}...\n")
    for seq in SeqIO.parse(seqsfile, 'fasta'):
        if seq.id==p.label.split(';')[1]:
            query=str(seq.seq).replace('*','')
            PA=ProteinAnalysis(query)
            gravy=PA.gravy()
            pc3top.at[ix,'Gravy']=gravy
    print('Done')
            
pc3top['hydrophobicity']=pc3top['Gravy'].apply(lambda row: 'hydrophobic' if row>0 else ('likely soluble' if row <0 else 'neutral'))

pc3top.to_csv('/'.join([wdir,'results/PCA_FITvsStool_PC3_loadings_top100.csv']),index=False)

# =============================================================================
# Plot values
# =============================================================================

def plot_interest(pc3top):
    
    pc3top=pc3top.sort_values(by='PC3',ascending=False).set_index('label')
    colors = ["#D85A30" if v > 0 else ("#378ADD" if v<0 else "#737373") for v in pc3top['Gravy'].tolist()]
    #red: hydrophobic; blue: likely soluble; grey: neutral
    
    s=pc3top['PC3']
    fig, ax = plt.subplots(figsize=(100 * 0.35 + 1,6))
    ax.vlines(s.index, 0, s.values, colors=colors, linewidth=2)
    ax.scatter(s.index, s.values,  color=colors, zorder=3, s=30)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel("PC3loading")
    ax.set_title("Top 100 protein loadings in PC3")
    ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
    sb.despine()
    plt.tight_layout()
    plt.show()
    
plot_interest(pc3top)
plt.savefig('/'.join([wdir,'results/Hydrophobicity_top100.pdf']))
