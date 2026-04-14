#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 16 10:21:35 2026

Spin protocol reproducibility

@author: ekateria
"""

import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator as annot
from itertools import combinations as comb

wdir='PATH_TO_YOUR_WORKDIR'

Data=pd.read_excel('/'.join([wdir,'results/combined_AllData_FITvsStool_bactonly.xlsx']))
Samples=pd.read_excel('/'.join([wdir,'results/combined_protein_summary_FITvsStool.xlsx']))
dbannot=pd.read_csv('/'.join([wdir.replace('BSA_removal',''), 'database/eggnog_annot.emapper.annotations']), sep='\t',skiprows=4)

# =============================================================================
# Calculate reproducibility (Fraction of consistently recovered proteins)
# =============================================================================

def reproducibility(data,pids,remtypes,thr, mode,Samples=Samples):
    data=data.copy()
    Repr=pd.DataFrame()
    for pid in pids:
        for remtype in remtypes:
            samnames=Samples.loc[Samples['BSA_removal']==remtype]
            samnames=samnames.loc[samnames['Pid']==pid,'SampleID'].unique().tolist()
            #remove BC samples
            samnames=[s for s in samnames if 'BC' not in s]
            samnames.append('Protein')
            d=data[samnames].set_index('Protein')
            d[d>thr]=1
            d['Sum']=d.sum(axis=1)
            d=d.query('Sum!=0')
            
            if mode=='all':
                numall=len(d.query(f'Sum>={len(samnames)-1}')) #in all replicates
            else:
                numall=len(d.query('Sum>=2')) #in all replicates

            
            reprod=pd.DataFrame({'Pid':[pid],'BSA_removal':[remtype],'NumSamples':[len(samnames)-1],
                                 'Total_NumProtDet':[len(d)],
                                 'Constistent_NumProtDet':[numall],
                                 'Reproducibility, %':[numall/len(d)*100]})
            Repr=pd.concat([Repr,reprod])
    
    return Repr
        
pids=['A','B','C']   
remtypes=['Spin','Stool_noBSA']
thr=0
mode='at least two' #mode='all'; count number of proteins detected in at least two/all replicates
Repr=reproducibility(Data,pids,remtypes, thr, mode)
Repr.to_excel('/'.join([wdir,f'results/Reproducibility_Spin_vs_Stool_thr{thr}_{mode}.xlsx']), index=False)

ax=sb.barplot(Repr,x='Reproducibility, %',y='BSA_removal',hue='Pid',palette='husl')
ax.set(xlabel=f"Fraction of proteins detected in {mode} replicates, %",ylabel='')
ax.set_yticklabels(['FIT','Stool'])

plt.savefig('/'.join([wdir,f'results/Reproducibility_by_replicates_{mode}.pdf']))

# =============================================================================
# Create venn diagrams for each individual
# =============================================================================

def reproducibility_by_pid(data, pids, thr, Samples=Samples):
    
    data=data.copy()   
    
    #make a proteins list
    def list_proteins(df,sam):
        
        protlist=[]
        for s in sam:
            d=df.loc[df[s]>0]
            pl=d.index.tolist()
            protlist.extend(pl)
            
        protlist=list(set(protlist))
        
        return protlist
    
    Venn_diagram=pd.DataFrame()
    
    Stoolonly={}
    FITonly={}
    Common={}
    for pid in pids:
        samnames=Samples.loc[Samples['Pid']==pid,'SampleID'].unique().tolist()
        #remove BC samples
        samnames=[s for s in samnames if 'BC' not in s]
        d=data[samnames]
        
        d[d<=thr]=0 #set all values below threshold to zero
        d['Sum']=d.sum(axis=1)
        d=d.query('Sum>0') #remove empty rows
        
        stsam=[s for s in samnames if 'stool' in s]
        fitsam=[s for s in samnames if 'stool' not in s]
        
        stlist=list_proteins(d,stsam)
        fitlist=list_proteins(d,fitsam)
        comlist=list(set(stlist).intersection(set(fitlist)))
        
        Stoolonly[pid]=[p for p in stlist if p not in comlist]
        FITonly[pid]=[p for p in fitlist if p not in comlist]
        Common[pid]=comlist
        
        vend=pd.DataFrame({'Pid':[pid],'TotalNumProt':[len(d)],'InStool':[len(Stoolonly[pid])],
                           'InFIT':[len(FITonly[pid])],'InCommon':[len(Common[pid])]})
        
        Venn_diagram=pd.concat([Venn_diagram,vend])
                          
    return Venn_diagram, Stoolonly, FITonly, Common

Data=Data.set_index('Protein')
VennDiag, Stoolonly, FITonly, Common =reproducibility_by_pid(Data.loc[Data.index.str.contains('Gene')],pids,thr)

VennDiag.to_excel('/'.join([wdir,f'results/VennDiag_Spin_vs_Stool_thr{thr}_MaxLFQ.xlsx']), index=False)

# =============================================================================
# Find relative abundance of those proteins that are shared/unique to stool
# and FIT samples
# =============================================================================

def find_prot_conc(data, colname, Stoolonly=Stoolonly, FITonly=FITonly, Common=Common, Samples=Samples):
    
    data=data.copy()
    
    AllData=pd.DataFrame()
    
    for pid in pids:
        samnames=Samples.loc[Samples['Pid']==pid,'SampleID'].unique().tolist()
        #remove BC samples
        samnames=[s for s in samnames if 'BC' not in s]
        d=data[samnames]
        d[d<=thr]=0 #set all values below threshold to zero
        
        st=d.loc[Stoolonly[pid]]
        st['Type']='Stoolonly'
        fit=d.loc[FITonly[pid]]
        fit['Type']='FITonly'
        common=d.loc[Common[pid]]
        common['Type']='Common'
        
        alldata=st.copy()
        alldata=pd.concat([alldata,fit])
        alldata=pd.concat([alldata,common])
        alldata['Pid']=pid
        alldata=alldata.reset_index()
        
        AllData=pd.concat([AllData,alldata])
    
    AllData=pd.melt(AllData, id_vars=['Protein','Type','Pid'], value_name=colname, var_name='Sample')
    AllData=AllData.dropna(subset=colname)
    
    return AllData
        
FinData=find_prot_conc(Data.loc[Data.index.str.contains('Gene')],'MaxLFQ')


hue_order=['FITonly','Common','Stoolonly']

pairs = []

for p in pids:
    for c1, c2 in comb(hue_order, 2):
        pairs.append(((p, c1), (p, c2)))

print(pairs)

ax=sb.boxplot(FinData.query('MaxLFQ>0'),x='Pid',y='MaxLFQ',hue='Type',
                 log_scale=True, order=pids, hue_order=hue_order, legend=False,
                 palette=['#ec9e7d','#9697b1','#7ebfaa'])

ann=annot(ax,pairs,data=FinData.query('MaxLFQ>0'),x='Pid',y='MaxLFQ', hue='Type', order=pids, hue_oder=hue_order)
ann.configure(test='Mann-Whitney',text_format='star',loc='inside',hide_non_significant=True)
ann.apply_and_annotate()   
ax.set(ylabel='MaxLFQ',xlabel='')
plt.savefig('/'.join([wdir,'results/Reproducibility_by_ID_ProtAbundance_MaxLFQ.pdf']))
        
# =============================================================================
# Find annotations of the proteins detected only in stool/FIT or in common
# =============================================================================

def format_genelists(qdict,descr):
    Prot=pd.DataFrame()
    for p in qdict:
        prot=pd.DataFrame(qdict[p],columns=['Protein'])
        prot['Pid']=p
        prot['DetectedIn']=descr
        Prot=pd.concat([Prot,prot])
        
    return Prot

fit=format_genelists(FITonly,'FITonly')
stool=format_genelists(Stoolonly,'Stoolonly')
common=format_genelists(Common,'Common')

DetSummary=pd.concat([fit, stool])
DetSummary=pd.concat([DetSummary,common])

dbannot=dbannot.rename(columns={'#query':'Protein'})

DetSummary=DetSummary.merge(dbannot[['Protein','COG_category','Preferred_name','Description','max_annot_lvl']])

#Make a summary table
# count proteins per (id, detection, COG_category)
counts = (DetSummary.groupby(['Pid', 'DetectedIn', 'COG_category']).size().reset_index(name='count'))
#Calculate fraction of each COG category relative to within each detection category (i.e. how many of stool only belong to a given COG)
totals = DetSummary.groupby(['Pid','DetectedIn']).size().reset_index(name='total_bydet')
result = counts.merge(totals, on=['Pid','DetectedIn'])
result['PercentByDet'] = result['count']/result['total_bydet']*100

# =============================================================================
# Make a plot - fraction of COGs in onyFIT/onlyStool/Common
# =============================================================================
Data=Data.reset_index()
cat_count=Data[['Protein']].merge(dbannot[['Protein','COG_category']])
cat_count=cat_count['COG_category'].value_counts().reset_index(name='ProtCount')

result=result.merge(cat_count,on='COG_category')
plot_res=result.query('count>10').sort_values(by='ProtCount',ascending=False)
plot_res=plot_res.loc[plot_res['COG_category']!='-']

ax=sb.boxplot(plot_res,hue='DetectedIn',hue_order=hue_order, y='PercentByDet',x='COG_category',
              legend=False, palette=['#ec9e7d','#9697b1','#7ebfaa'])
labels=pd.DataFrame([tick.get_text() for tick in ax.get_xticklabels()], columns=['COG_category'])
labels=labels.merge(cat_count, on='COG_category')
labels['label']=labels.apply(lambda row: f'{row.COG_category}; (Nprot={row.ProtCount})', axis=1)
ax.set_xticklabels(labels['label'].tolist(),rotation=90)
ax.set(ylabel='Proportion of proteins within COG category,%',xlabel='')

#Make separation lines between categories
n = len(ax.get_xticks())
for x in range(n - 1):
    ax.axvline(x + 0.5, linestyle='--', color='gray', linewidth=0.8, alpha=0.7)
    
plt.savefig('/'.join([wdir,'results/Fraction_by_COG_and_DetectedIn_hor.pdf']))

pairs = []

for c in labels['COG_category'].tolist():
    print(c)
    pl=plot_res.loc[plot_res['COG_category']==c]
    if len(pl)>3:
        if len(pl['Pid'].unique().tolist())>=2 and len(pl['DetectedIn'].unique().tolist())>1:
            comblist=pl['DetectedIn'].unique().tolist()
            for c1,c2 in comb(comblist,2):
                pairs.append(((c, c1), (c, c2)))

print(pairs)

ann=annot(ax,pairs,data=plot_res,x='COG_category',y='PercentByDet', hue='DetectedIn', hue_oder=hue_order)
ann.configure(test='Mann-Whitney',text_format='star',loc='inside',hide_non_significant=True)

ax, annotations = ann.apply_and_annotate()

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

MannWStats.to_csv('/'.join([wdir,'results/Fraction_by_COG_and_DetectedIn_MannW.csv']),sep='\t', index=False)


