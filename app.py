#from Bio import Entrez
#import ssl
#import urllib.request
import streamlit as st
import pandas as pd
import numpy as np
import re
st.title("Get flat file, sample file, and entrez list for lymphgen")
fileuse=st.file_uploader("Upload your file", type="txt",accept_multiple_files=True)
vaf = st.number_input("Enter VAF cutoff", min_value=0, max_value=100, value=8, step=1)
# Create an SSL context that doesn't verify the certificates

#ssl._create_default_https_context = ssl._create_unverified_context


@st.cache_data
def convert_df_to_csv(df):
    return df.to_csv(index=False,sep='\t').encode('utf-8')


def get_entrez_id(pd_gene_names):
    gene_ids=pd.read_csv("names_to_entrezids.txt", header=0,sep='\t',index_col=None)

    pd_gene_names=pd.merge(pd_gene_names,gene_ids,how='inner',left_on='GENE',right_on='external_gene_name')
    return pd_gene_names
if fileuse != []:
    st.write("File(s) uploaded")
    st.write("VAF cutoff is "+str(vaf))
    all_fs=[]
    for i in fileuse:
        f_use=pd.read_csv(i, sep="\t", header=0, index_col=0)
        
        f_use.reset_index(inplace=True)
        #breakpoint()
        for index, row in f_use.iterrows():
            if row['AF'] < vaf:
                f_use.drop(index, inplace=True)
            #if row['p.CHANGE'] == 'p.?':
            #    f_use.drop(index, inplace=True)
            #elif row['p.CHANGE'] == '':
            #    f_use.drop(index, inplace=True)
            elif row['EFFECT'] == 'structural_interaction_variant':
                f_use.drop(index, inplace=True)   
            elif re.findall('rs',str(row['rsID'])) and not re.findall('CO',str(row['COSMIC'])):
                f_use.drop(index, inplace=True)                
            #elif row['VARIANT_TYPE'] == 'SNP':
            #    f_use.drop(index, inplace=True)
            elif row['FILTER'] != 'PASS':   
                f_use.drop(index, inplace=True) 
        f_use=get_entrez_id(f_use)
        f_use['Sample']=f_use['#SAMPLE_NAME'].astype(str)
        f_use['ENTREZ.ID']=f_use['entrezgene_id']
        f_use['Type']=''
        f_use['Type'][f_use['EFFECT'].str.contains('synon')]='Synon'
        f_use['Type'][f_use['EFFECT'].str.contains('frameshift')]='TRUNC'
        f_use['Type'][f_use['EFFECT'].str.contains('stop_gain')]='TRUNC'
        f_use['Type'][f_use['EFFECT'].str.contains('splice_donor')]='TRUNC'
        f_use['Type'][f_use['EFFECT'].str.contains('splice_acceptor')]='TRUNC'
        f_use['Type'][f_use['Type'].str.contains('splice_region')]='TRUNC'
        f_use['Type'][f_use['EFFECT'].str.contains('missense')]='MUTATION'
        f_use['Type'][f_use['EFFECT'].str.contains('inframe')]='MUTATION'
        f_use['Type'][f_use['EFFECT'].str.contains('in_frame')]='MUTATION'
        if f_use['p.CHANGE'].str.contains('L265P').any():
            f_use['Type'][f_use['p.CHANGE'].str.contains('L265P')]='L265P'
        f_use_for_lymphgen=f_use[f_use['Type'] !='']
        all_fs.append(f_use_for_lymphgen)
    all_fs=pd.concat(all_fs,ignore_index=True)
    fordl = convert_df_to_csv(all_fs[['Sample','ENTREZ.ID','Type']])
    entrezid = convert_df_to_csv(pd.read_csv('entrezids.txt',sep='\t', header=0, index_col=None))
    sample_names=all_fs[['Sample']].drop_duplicates()
    st.write(all_fs)
    annotation_samples=convert_df_to_csv(pd.DataFrame({'Sample.ID':sample_names['Sample'],'Copy.Number':0, 'BCL2.transloc':'NA','BCL6.transloc':'NA'}))
    st.download_button("Download flat mutation file for lymphgen", fordl, "flatfilewmutations.txt", "text/csv", key='download-txt-mutants')
    st.download_button("Download entrez ids for lymphgen", entrezid, "entrezids.txt", "text/csv", key='download-txt-entrez')
    st.download_button("Download annotation file for lymphgen", annotation_samples, "sample_annotation.txt", "text/csv", key='download-txt-samples')    
