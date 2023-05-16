import pandas as pd
import math
import os
import numpy as np
import argparse

'''
command line
python3 get_gatk_cnv_long_format.py --list wxs_male_vs_male_delme.gatk_cnv.funcotated.tsv.gene_list.txt 
                --seg wxs_male_vs_male_delme.tumor.gatk_cnv.called.seg 
                --absolute_CN 2 
                --ploidy 4
'''

# Instantiate the parser
parser = argparse.ArgumentParser(description='Get continuous & discrete cnv')

# Required positional argument
parser.add_argument('--list', help='List of cnv file ')
parser.add_argument('--seg', help='Seg file for cnv ')
parser.add_argument('--absolute_CN',type=int, help='Seg file for cnv ')
parser.add_argument('--ploidy', type= int,help='Seg file for cnv ')
args = parser.parse_args()

absolute_CN=args.absolute_CN
ploidy=args.ploidy
cnv_dataframe= pd.read_csv(args.list,sep='\t',usecols = ['gene','segment_contig','segment_start','segment_end','segment_mean'], low_memory = True)
cnv_dataframe['cnv_score_continuous']=round(pow(2,cnv_dataframe['segment_mean']))

def define_score_metrix(absolute_CN,ploidy): 
    # Converting score to GITSIC style values
    #
    score_matrix={0:-absolute_CN}
    for i in range(1,absolute_CN+ploidy+1,1):
        if i < absolute_CN:
            score_matrix[i]=-1
        elif i ==absolute_CN:
            score_matrix[i]=0    
        else:
            score_matrix[i]=1   
    return score_matrix        

def convert(x):
    score_metrix=define_score_metrix(absolute_CN,ploidy)
    if x in score_metrix:
        x=score_metrix[x]
    else:
        x= 2
    return x    

def get_set_data():
    seg_data=[]
    with open(args.seg, "r") as f:
        switch_on="OFF"
        for line in f:
            #print (line)
            if "@RG" in line:
                line=line.replace('\n','')
                line_data=line.split("\t")
                sample_id=line_data[2].split(":")[1]
                #print(sample_id)
            if "START" in line or switch_on =="ON":
                switch_on ="ON"
                line=line.replace('\n','')
                line_data=line.split("\t")
                seg_data.append(line_data)
    return seg_data,sample_id

def map_genes_seg_to_ref(seg_file_data):
    header=seg_file_data[0]
    seg_file_data.pop(0)
    ref_genes=[]
    for calls in seg_file_data:
        gene=''
        for ref in ref_gene_bed_list:
            if (
                ref[0] == calls[0] and  #chr
                int(ref[1]) < int(calls[1]) and  #start
                int(ref[2]) > int(calls[2])
                ): #end
                gene=ref[3]
                break
        ref_genes.append(gene) 

    seg_file_data_df=pd.DataFrame(seg_file_data,columns=header)
    seg_file_data_df[header[1]]=seg_file_data_df[header[1]].apply(int)
    seg_file_data_df[header[2]]=seg_file_data_df[header[2]].apply(int)
    seg_file_data_df[header[3]]=seg_file_data_df[header[3]].apply(int)
    seg_file_data_df[header[4]]=seg_file_data_df[header[4]].apply(float)
    seg_file_data_df['gene']=ref_genes
         
    return seg_file_data_df  

cnv_dataframe["cnv_score_discrete"]=cnv_dataframe["cnv_score_continuous"].map(convert)

#path to ref file to map genes
path_current_file=str(os.path.join(os.path.dirname(__file__))).split("/")
path_current_file[-1]='REFS'
path_current_file.append('Homo_sapiens.GRCh38.105.chr.gtf_genes.bed')
path_ref_file='/'.join(path_current_file)

ref_gene_bed=pd.read_csv(path_ref_file,sep='\t',names=["chr",'start','end','gene'])
ref_gene_bed['chr']='chr'+ref_gene_bed['chr']
ref_gene_bed_list=ref_gene_bed.values.tolist()

seg_data,sample_ID=get_set_data()
seg_data_df=map_genes_seg_to_ref(seg_data)
#seg_data_df.to_csv("test",sep='\t')
seg_data_df.sort_values(by=['CONTIG','START'],inplace=True)
seg_data_df=seg_data_df.convert_dtypes()

cnv_dataframe.drop(['gene'],axis=1,inplace=True)
cnv_dataframe.drop_duplicates(inplace=True)
cnv_dataframe.sort_values(by=['segment_contig','segment_start'],inplace=True)
merged_seg_cnv=pd.merge(cnv_dataframe,seg_data_df,how='inner',right_on='START', left_on='segment_start')
merged_seg_cnv.drop(["segment_start","segment_end","CONTIG"],axis=1,inplace=True)
merged_seg_cnv.rename(columns={"segment_contig":"chr","gene":"Hugo_Symbol"},inplace=True)
merged_seg_cnv["SAMPLE_ID"]=sample_ID
merged_seg_cnv=merged_seg_cnv[["SAMPLE_ID","Hugo_Symbol","chr","START","END","segment_mean","cnv_score_continuous","cnv_score_discrete","NUM_POINTS_COPY_RATIO","MEAN_LOG2_COPY_RATIO"]]
#print(merged_seg_cnv)
merged_seg_cnv.to_csv("per_sample.txt",sep='\t',index=False)

cnv_discrete_per_sample = merged_seg_cnv[["Hugo_Symbol","chr","START","END","cnv_score_discrete"]]
cnv_discrete_per_sample=cnv_discrete_per_sample.rename(columns={"cnv_score_discrete":str(sample_ID)})
#cnv_discrete_per_sample['Entrez_Gene_Id']=''
cnv_discrete_per_sample=cnv_discrete_per_sample[['Hugo_Symbol',str(sample_ID)]]
#cnv_discrete_per_sample['Hugo_Symbol'].replace('', np.nan,inplace=True)
cnv_discrete_per_sample=cnv_discrete_per_sample[cnv_discrete_per_sample["Hugo_Symbol"].astype(bool)]
#cnv_discrete_per_sample.dropna(subset=['Hugo_Symbol'], inplace=True)
cnv_discrete_per_sample.to_csv("cnv_discrete_per_sample.txt",sep='\t',index=False)

cnv_continuous_per_sample = merged_seg_cnv[["Hugo_Symbol","chr","START","END","cnv_score_continuous"]]
cnv_continuous_per_sample=cnv_continuous_per_sample.rename(columns={"cnv_score_continuous":str(sample_ID)})
cnv_continuous_per_sample=cnv_continuous_per_sample[cnv_continuous_per_sample["Hugo_Symbol"].astype(bool)]
#cnv_continuous_per_sample['Entrez_Gene_Id']=''
cnv_continuous_per_sample=cnv_continuous_per_sample[['Hugo_Symbol',str(sample_ID)]]
cnv_continuous_per_sample.to_csv("cnv_continuous_per_sample.txt",sep='\t',index=False)


#cnv_dataframe_gene=cnv_dataframe[['gene','segment_contig','segment_start','segment_end']].values.tolist()
#print(len(cnv_dataframe_gene))

#chunks=16
#list_df = [
#        cnv_dataframe[i : i + chunks].copy(deep=True) for i in range(0, len(cnv_dataframe), chunks)
#    ]
#print(len(list_df))
#genes=map_genes(cnv_dataframe_gene)

'''
def pick_gene(row):
    print(row["segment_start"] > ref_gene_bed["start"] and ref_gene_bed['chr'] == row["segment_contig"])
    return 1#[ref_gene_bed['chr'] == row["segment_contig"] & int(row["segment_start"]) > int(ref_gene_bed["start"])]
'''
'''
def searching(x):
    ap=[]
    for start,end,gene in zip(ref_gene_bed['start'],ref_gene_bed['end'],ref_gene_bed['gene']):
        if start < x:# and x < end :
            if x > end:
                ap.append(True)
        else:
            ap.append(False)
    return ap 
'''     
'''
def map_genes(cnv_data):
    ref_genes=[]
    cnv_data_list=cnv_data.values.tolist()
    for calls in cnv_data_list:
        gene=''
        for ref in ref_gene_bed_list:
            if ref[0] == calls[1] and ref[1] < calls[2] and ref [2] > calls[3]:
                gene=ref[3]
                break
        ref_genes.append(gene) 

    #cnv_data.drop(columns=['gene'],axis=1,inplace=True)
    cnv_data['gene1']=ref_genes
         
    return cnv_data    
'''
'''
from concurrent.futures import ThreadPoolExecutor
threads=[]
results=[]
with ThreadPoolExecutor(max_workers=10) as executor:
    for i in list_df:
        threads.append(executor.submit(map_genes,i))

    for t in threads:
        results.append(t.result())   
df = pd.concat(results)         
#print(df)
'''
'''
for calls in cnv_dataframe_gene:
    gene=''
    for ref in ref_gene_bed_list:
        if ref[0] == calls[1] and ref[1] < calls[2] and ref [2] > calls[3]:
            gene=ref[3]
            break
    ref_genes.append(gene)    
'''    