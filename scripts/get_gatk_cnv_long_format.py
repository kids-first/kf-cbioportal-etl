import pandas as pd
import math
import os
import numpy as np
import argparse
from concurrent.futures import ThreadPoolExecutor
import sys
'''
command line
python3 get_gatk_cnv_long_format.py --segfiles_folder_path gatk_copy_ratio_segments_tumor --manifest manifest.tsv

'''

# Instantiate the parser
parser = argparse.ArgumentParser(description='Get continuous & discrete cnv')

# Required positional argument
parser.add_argument('--segfiles_folder_path', help='provide seg file folder to read seg files')
parser.add_argument('--manifest', help='provide manifest file to map BS IDs with file output in tsv format')
#parser.add_argument('--seg', help='Seg file for cnv ')
#parser.add_argument('--absolute_CN',type=int, help='Seg file for cnv ')
#parser.add_argument('--ploidy', type= int,help='Seg file for cnv ')
args = parser.parse_args()

absolute_CN=2#args.absolute_CN
ploidy=4#args.ploidy
#cnv_dataframe= pd.read_csv(args.list,sep='\t',usecols = ['gene','segment_contig','segment_start','segment_end','segment_mean'], low_memory = True)
#cnv_dataframe['cnv_score_continuous']=round(pow(2,cnv_dataframe['segment_mean']))

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

def get_seg_data(seg):
    seg_data=[]
    with open(seg, "r") as f:
        switch_on="OFF"
        for line in f:
            #print (line)
        #    if "@RG" in line:
        #        line=line.replace('\n','')
        #        line_data=line.split("\t")
        #        sample_id=line_data[2].split(":")[1]
                #print(sample_id)
            if "START" in line or switch_on =="ON":
                switch_on ="ON"
                line=line.replace('\n','')
                line_data=line.split("\t")
                seg_data.append(line_data)
    return seg_data

def map_genes_seg_to_ref(seg_file_data,ref_gene_bed_list):
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

def get_ref_bed_data(ref_file_name):

    path_current_file=str(os.path.join(os.path.dirname(__file__))).split("/")
    path_current_file[-1]='REFS'
    path_current_file.append('Homo_sapiens.GRCh38.105.chr.gtf_genes.bed')
    path_ref_file='/'.join(path_current_file)
    ref_gene_bed=pd.read_csv(path_ref_file,sep='\t',names=["chr",'start','end','gene'])
    ref_gene_bed['chr']='chr'+ref_gene_bed['chr']
    ref_gene_bed_list=ref_gene_bed.values.tolist()

    return ref_gene_bed_list

#path to ref file to map genes
def prepare_gat_cnv_from_seg_file(seg_file,sample_ID,ref_gene_bed_list):
    
    seg_data=get_seg_data(seg_file)
    seg_data_df=map_genes_seg_to_ref(seg_data,ref_gene_bed_list)
    seg_data_df.sort_values(by=['CONTIG','START'],inplace=True)
    seg_data_df=seg_data_df.convert_dtypes()
    seg_data_df['cnv_score_continuous']=round(pow(2,seg_data_df['MEAN_LOG2_COPY_RATIO']))
    seg_data_df["cnv_score_discrete"]=seg_data_df["cnv_score_continuous"].map(convert)
    seg_data_df=seg_data_df.drop(["MEAN_LOG2_COPY_RATIO","CALL","NUM_POINTS_COPY_RATIO",'CONTIG','START','END'],axis=1)
    seg_data_df=seg_data_df[seg_data_df["gene"].astype(bool)]
    cnv_per_sample=seg_data_df
    cnv_per_sample.rename(columns={"gene":"Hugo_Symbol"},inplace=True)

    df_without_duplicate_gene = cnv_per_sample[cnv_per_sample['Hugo_Symbol'].duplicated(keep=False) == False] # dropping all the duplicate genes
    df_duplicate_genes = cnv_per_sample[cnv_per_sample['Hugo_Symbol'].duplicated(keep=False) == True] # extracting duplicate genes
    df2_loss_cnv=df_duplicate_genes[df_duplicate_genes['cnv_score_discrete']<0] #creteria for loss
    df2_gain_cnv=df_duplicate_genes[df_duplicate_genes['cnv_score_discrete']>0] #creteria for gain
    df2_loss_cnv = df2_loss_cnv.sort_values('cnv_score_discrete', ascending=True).drop_duplicates('Hugo_Symbol').sort_index() # remove duplicate gene and pick one with min cnv score for loss
    df2_gain_cnv = df2_gain_cnv.sort_values('cnv_score_discrete', ascending=False).drop_duplicates('Hugo_Symbol').sort_index() # remove duplicate gene and pick one with max cnv score for gain
    
    df_without_duplicate_gene=pd.concat([df_without_duplicate_gene,df2_loss_cnv,df2_gain_cnv]) #Add duplicate gene back
    #print(df_without_duplicate_gene)
    #print(df2_loss_cnv)
    #print(df2_gain_cnv)
    #print(df_without_duplicate_gene)
    cnv_discrete_per_sample=df_without_duplicate_gene.drop(["cnv_score_continuous"],axis=1) #prepare df for discrete scores
    cnv_discrete_per_sample=cnv_discrete_per_sample.rename(columns={"cnv_score_discrete":str(sample_ID)})
    cnv_continuous_per_sample=df_without_duplicate_gene.drop(['cnv_score_discrete'],axis=1) #prepare df for continuous scores
    cnv_continuous_per_sample['cnv_score_continuous']=cnv_continuous_per_sample["cnv_score_continuous"].astype('int')
    cnv_continuous_per_sample=cnv_continuous_per_sample.rename(columns={'cnv_score_continuous':str(sample_ID)})

    return [cnv_discrete_per_sample,cnv_continuous_per_sample]

def outer_merge_list(df_list):
    merged_list=df_list[0]
    for i in range(1,len(df_list),1):
        merged_list=pd.merge(merged_list,df_list[i],how='outer',on=["Hugo_Symbol"])
    return merged_list

if not os.path.isfile(args.manifest):
        raise TypeError("Cannot find given manifest file") 
else:
    manifest=pd.read_csv(args.manifest,sep='\t',header=0)
    file_name=manifest['s3_path'].str.split("/",expand=True)[5]
    file_name=file_name.rename("filename")
    map_bs_filename=pd.concat([manifest['biospecimen_id'],file_name],axis=1)

ref_gene_bed_list=get_ref_bed_data('Homo_sapiens.GRCh38.105.chr.gtf_genes.bed')   #used for gene mapping to seg calls 

threads=[]
results=[]
counter=0
with ThreadPoolExecutor(max_workers=100) as executor:
    for filename in os.listdir(args.segfiles_folder_path):
        file_path = os.path.join(args.segfiles_folder_path, filename)
        # checking if it is a file
        if os.path.isfile(file_path) and file_path.endswith(".seg"):
            file_name=file_path.split('/')[1] # file name from split
            sample_ID=map_bs_filename[map_bs_filename["filename"]==file_name]['biospecimen_id'].values[0]
            file_path='wxs_male_vs_male_delme.tumor.gatk_cnv.called.seg'
            threads.append(executor.submit(prepare_gat_cnv_from_seg_file,file_path,sample_ID,ref_gene_bed_list))
            counter=counter+1
            if counter == 10:
                break
    for active_threads in threads:
        results.append(active_threads.result())              

# code to split it into 2 lists
discrete_list, continous_list = map(list, zip(*results))
df_discrete=outer_merge_list(discrete_list)
df_discrete=df_discrete.replace(np.nan,'', regex=True)
df_continous=outer_merge_list(continous_list)
df_continous=df_continous.replace(np.nan,'', regex=True)
df_continous=df_continous.sort_values("Hugo_Symbol")
df_discrete=df_discrete.sort_values("Hugo_Symbol")
df_continous.to_csv("cnv_continous.txt",sep='\t',index=False)
df_discrete.to_csv("cnv_discrete.txt",sep='\t',index=False)
#print(df_discrete)
#print(df_discrete)
#print(df_continous)

#print(cnv_discrete_per_sample)
#print(cnv_continuous_per_sample)

'''
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
'''

#cnv_dataframe['cnv_score_continuous']=round(pow(2,cnv_dataframe['segment_mean']))
#cnv_dataframe["cnv_score_discrete"]=cnv_dataframe["cnv_score_continuous"].map(convert)

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