import pandas as pd
import math
import os
import numpy as np
import argparse
from concurrent.futures import ThreadPoolExecutor
import sys
import json

'''
command line
python3 get_gatk_cnv_long_format.py --segfiles_folder_path gatk_copy_ratio_segments_tumor --manifest manifest.tsv

'''

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
    absolute_CN = 2 
    ploidy = 4
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
    #read ref data to map genes

    path_current_file=str(os.path.join(os.path.dirname(__file__))).split("/")
    path_current_file[-1]='REFS'
    path_current_file.append(gene_bed_ref_filename)
    path_ref_file='/'.join(path_current_file)
    ref_gene_bed=pd.read_csv(path_ref_file,sep='\t',names=["chr",'start','end','gene'])
    ref_gene_bed['chr']='chr'+ref_gene_bed['chr']
    ref_gene_bed_list=ref_gene_bed.values.tolist()

    return ref_gene_bed_list

def extract_integer(string_data):
    return string_data.replace('chr','')

def prepare_seg_output(seg_data,cbio_ID):
    #seg_data_df=seg_data_df.drop(["CALL"])
    seg_data["CONTIG"]=seg_data["CONTIG"].map(extract_integer)
    seg_data=seg_data.drop(['CALL','gene'],axis=1)
    seg_data["ID"]=str(cbio_ID)
    seg_data=seg_data.rename(columns={"CONTIG":"chrom","START":"loc.start","END":"loc.end","NUM_POINTS_COPY_RATIO":"num.mark","MEAN_LOG2_COPY_RATIO":"seg.mean"})
    ID_i = seg_data.pop("ID")
    seg_data.insert(0, 'ID', ID_i)
    return seg_data


def prepare_gat_cnv_from_seg_file(seg_file,sample_ID,ref_gene_bed_list,map_cbio_BS_IDs):
    #prepare data in wide format for per seg file

    seg_data=get_seg_data(seg_file)
    seg_data_df=map_genes_seg_to_ref(seg_data,ref_gene_bed_list)
    seg_data_df.sort_values(by=['CONTIG','START'],inplace=True)
    seg_data_df=seg_data_df.convert_dtypes()

    #print(sample_ID)
    #print(seg_file)
    cbio_ID_list=map_cbio_BS_IDs[map_cbio_BS_IDs['T_CL_BS_ID']==sample_ID]['Cbio_Tumor_Name'].values
    if len(cbio_ID_list)>0:
        cbio_ID=cbio_ID_list[0]
    else:    
        print("cbio tumor id not found for ", sample_ID, "therefore ignoring this sample")
        return []
        
    #cbio_ID=sample_ID

    seg_processing=seg_data_df.copy()
    seg_output=prepare_seg_output(seg_processing,cbio_ID)
    
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

    #print()
    
    cnv_discrete_per_sample=df_without_duplicate_gene.drop(["cnv_score_continuous"],axis=1) #prepare df for discrete scores
    cnv_discrete_per_sample=cnv_discrete_per_sample.rename(columns={"cnv_score_discrete":str(cbio_ID)})
    cnv_continuous_per_sample=df_without_duplicate_gene.drop(['cnv_score_discrete'],axis=1) #prepare df for continuous scores
    cnv_continuous_per_sample['cnv_score_continuous']=cnv_continuous_per_sample["cnv_score_continuous"].astype('int')
    cnv_continuous_per_sample=cnv_continuous_per_sample.rename(columns={'cnv_score_continuous':str(cbio_ID)})


    return [cnv_discrete_per_sample,cnv_continuous_per_sample,seg_output]

def outer_merge_list(df_list,on_header):
    merged_list=df_list[0]
    for i in range(1,len(df_list),1):
        merged_list=pd.merge(merged_list,df_list[i],how='outer',on=on_header)
    return merged_list

def run_thread_pool():
    threads=[]
    results=[]
    counter=0
    with ThreadPoolExecutor(max_workers=100) as executor:
        for filename in os.listdir(segfiles_folder_path):
            file_path = os.path.join(segfiles_folder_path, filename)
            # checking if it is a file
            if os.path.isfile(file_path) and file_path.endswith(".seg"):
                file_name=file_path.split('/')[-1] # extract file name from split
                sample_ID=map_bs_filename[map_bs_filename["filename"]==file_name]['biospecimen_id'].values[0]
                #file_path='wxs_male_vs_male_delme.tumor.gatk_cnv.called.seg'
        #        sample_ID='BS_0GS82Q3E'
        #        file_path='/Users/phuls/tools/kf-cbioportal-etl/scripts/gatk_copy_ratio_segments_tumor/ddc5c5ea-399c-4ded-9fce-48a275d9bbb4.tumor.gatk_cnv.called.seg'
                threads.append(executor.submit(prepare_gat_cnv_from_seg_file,file_path,sample_ID,ref_gene_bed_list,map_cbio_BS_IDs))
                
                counter=counter+1
                if counter == 100:
                    break
        for active_threads in threads:
            results.append(active_threads.result())        
        #        results.append(prepare_gat_cnv_from_seg_file(file_path,sample_ID,ref_gene_bed_list,map_cbio_BS_IDs))
        results_not_empty = [x for x in results if len(x) > 0 ]
    return results_not_empty  
import concurrent.futures

def run_process_pool():
    process=[]
    results=[]
    counter=0
    with concurrent.futures.ProcessPoolExecutor(max_workers=12) as executor:
        for filename in os.listdir(segfiles_folder_path):
            file_path = os.path.join(segfiles_folder_path, filename)
            # checking if it is a file
            if os.path.isfile(file_path) and file_path.endswith(".seg"):
                file_name=file_path.split('/')[-1] # extract file name from split
                sample_ID=map_bs_filename[map_bs_filename["filename"]==file_name]['biospecimen_id'].values[0]
                #file_path='wxs_male_vs_male_delme.tumor.gatk_cnv.called.seg'
        #        sample_ID='BS_0GS82Q3E'
        #        file_path='/Users/phuls/tools/kf-cbioportal-etl/scripts/gatk_copy_ratio_segments_tumor/ddc5c5ea-399c-4ded-9fce-48a275d9bbb4.tumor.gatk_cnv.called.seg'
                process.append(executor.submit(prepare_gat_cnv_from_seg_file,file_path,sample_ID,ref_gene_bed_list,map_cbio_BS_IDs))
                
               # counter=counter+1
               # if counter == 100:
               #     break
        for active_threads in process:
            results.append(active_threads.result())        
        #        results.append(prepare_gat_cnv_from_seg_file(file_path,sample_ID,ref_gene_bed_list,map_cbio_BS_IDs))
        results_not_empty_list = [x for x in results if len(x) > 0 ]
    return results_not_empty_list     

if __name__=="__main__":
        # Instantiate the parser
    parser = argparse.ArgumentParser(description='Get continuous & discrete cnv')

    # Required positional argument
    parser.add_argument('--config', help='Seg file for cnv ')
    args = parser.parse_args()

    # Opening JSON file
    json_file = open(args.config)
    json_dict = json.load(json_file)

    #read variables from json dict
    segfiles_folder_path=json_dict['file_loc_defs']['gatk_seg']['segfiles_folder_path']
    manifest_path=json_dict['file_loc_defs']['gatk_seg']['manifest']
    gene_bed_ref_filename=json_dict['file_loc_defs']['gatk_seg']['gene_bed_ref_filename']
    output_folder_path=json_dict['file_loc_defs']['gatk_seg']['output_folder_path']
    #absolute_CN=json_dict['file_loc_defs']['gatk_seg']['absolute_CN']
    #ploidy=json_dict['file_loc_defs']['gatk_seg']['ploidy_gain']

    if not os.path.isfile(manifest_path):
        raise TypeError("Cannot find given manifest file") 
    else:
        manifest=pd.read_csv(manifest_path,sep='\t',header=0)
        file_name=manifest['s3_path'].str.split("/",expand=True)[5]
        file_name=file_name.rename("filename")
        map_bs_filename=pd.concat([manifest['biospecimen_id'],file_name],axis=1)

    map_cbio_BS_IDs=pd.read_csv("pbta_all_genomics_etl_file.csv")
    map_cbio_BS_IDs=map_cbio_BS_IDs[map_cbio_BS_IDs["T_CL_BS_ID"].duplicated()==False]

    ref_gene_bed_list=get_ref_bed_data(gene_bed_ref_filename)   #used for gene mapping to seg calls 

    #results=run_thread_pool()
    results=run_process_pool()

    # code to split results list into 2 separeatelists
    discrete_list, continous_list,seg_list = map(list, zip(*results))
    df_discrete=outer_merge_list(discrete_list,"Hugo_Symbol")
    df_discrete=df_discrete.replace(np.nan,'', regex=True)
    df_continous=outer_merge_list(continous_list,"Hugo_Symbol")
    df_continous=df_continous.replace(np.nan,'', regex=True)
    df_continous=df_continous.sort_values("Hugo_Symbol")
    df_discrete=df_discrete.sort_values("Hugo_Symbol")

    df_seg=pd.concat(seg_list)
    df_seg=df_seg.replace(np.nan,'', regex=True)
    df_seg=df_seg.sort_values(["ID","chrom","loc.start"])

    df_continous.to_csv(os.path.join(output_folder_path,'cnv_continous.txt'),sep='\t',index=False)
    df_discrete.to_csv(os.path.join(output_folder_path,"cnv_discrete.txt"),sep='\t',index=False)
    df_seg.to_csv(os.path.join(output_folder_path,"cnv_seg.txt"),sep='\t',index=False)
    #print(df_seg)
    #print(df_discrete)
    #print(df_continous)
