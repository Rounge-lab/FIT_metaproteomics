# -*- coding: utf-8 -*-
"""

Preprocess metaproteomic tables (keep only proteins with probability >0.9 and >=2 peptides);
Calculate fractions of microbial/human/BSA proteins; 
Keep microbial proteins only

"""

import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt


wdir='PATH_TO_YOUR_WORKING_DIRECTORY'

def format_table(file,filtprot):
    data=pd.read_excel('/'.join([wdir,f'input/{file}']))
    
    dcols=data.columns.tolist()
    
    keep=dcols[0:11] #Keep columns with all the Protein characteristics
    
    samlist=[d for d in dcols if 'MaxLFQ Intensity' in d] #get the list of all samples from that run
    
    keep.extend(samlist) #add samples intensity
    keep.append(dcols[-1]) #add indistinguishable proteins
    
    data=data[keep]
    
    run=file.split('_')[-1].split('.')[0]
    newnames=[s.replace('MaxLFQ Intensity',run) for s in samlist] #rename the samples
    ids=[n[0] if 'NC' not in n else 'NC' for n in newnames] #get Pid names
    Run=[run for i in range(0,len(ids),1)]#add the run name

    samples=pd.DataFrame({'SampleID':newnames,'Run':Run, 'Pid':ids})
    
    data=data.rename(columns=dict(zip(samlist,newnames)))
    
    if filtprot:
        data=data.loc[data['Protein Probability']>0.9]
        data=data.loc[data['Combined Total Peptides']>=2]    
    
    return data,samples

def find_fractions(data,samids):

    bsa=data.loc[data['Protein'].str.contains('ALBU_BOVIN'),samids]/data[samids].sum(axis=0)*100
    bsa['type']='BSA'
    
    human=data.loc[data['Organism']=='Homo sapiens',samids].sum(axis=0)/data[samids].sum(axis=0)*100
    human=pd.DataFrame(human).T
    human['type']='human'

    bact=data.loc[data['Organism'].isna(),samids].sum(axis=0)/data[samids].sum(axis=0)*100
    bact=pd.DataFrame(bact).T
    bact['type']='bact'

    intensity=pd.DataFrame(data[samids].sum(axis=0)).T
    intensity.index=['total_intensity']

    summary=pd.concat([bsa,human])
    summary=pd.concat([summary,bact])
    summary=summary.set_index('type')

    summary=pd.concat([summary,intensity])
    
    return summary

# =============================================================================
# Make a summary file with total intensity and fractions
# =============================================================================

file='simulated_FragPipe_output.xlsx'
Data,Samples=format_table(file,True)
Summary=find_fractions(Data,Samples['SampleID'].tolist())

Samples['Type']=Samples['SampleID'].apply(lambda row: 'Stool' if 'stool' in row else ('NC' if 'N' in row else ('FIT' if 'BC' not in row else 'Bead-beating')))
Summary=Summary.T.reset_index().rename(columns={'index':'SampleID'})
Samples=Samples.merge(Summary,on='SampleID',how='left')
Samples=Samples.melt(id_vars=['SampleID','Run','Pid','Type','total_intensity'],var_name='ProteinType',value_name='Percentage')

Samples['BSA_removal']=Samples['Type'].apply(lambda row: 'Spin' if row=='FIT' else
                                             ('Spin+Beadbeating' if row=='Bead-beating' else ('Stool_noBSA' if row=='Stool' else ('Unprocessed' if row=='no BSA rem' else row))))

Samples['Type']=Samples['SampleID'].apply(lambda row: 'NC' if 'N' in row else ('Stool' if 'stool' in row else 'FIT'))
Samples.loc[Samples['SampleID']=='NC','Type']='FITbuffer'
Samples.loc[Samples['Pid'].str.contains('P'),'BSA_removal']='Unprocessed'
Samples.to_excel('/'.join([wdir,'results/combined_protein_summary.xlsx']), index=False)

# =============================================================================
# Plot fraction data
# =============================================================================
order=['FITbuffer','Unprocessed','Spin','Stool_noBSA','NC']
ax=sb.boxplot(Samples.loc[Samples['Type']!='NC'], x='BSA_removal',y='Percentage', hue='ProteinType', 
           order=order, palette=['#0097b2','#cca300','#cc4e00'],
           hue_order=['bact','human','BSA'], legend=False)
ax.set(xlabel='BSA removal strategy', ylabel='Fraction of Max_LFQ intensity from proteins,%')

plt.savefig(('/'.join([wdir,'results/bact_human_BSA_percentage.pdf'])))

AllData_bactonly=Data.loc[Data['Protein'].str.contains('Gene')]
AllData_bactonly.to_excel('/'.join([wdir,'input/true_MaxLFQ_nohuman.xlsx']), index=False)

 