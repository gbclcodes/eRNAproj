import numpy as np
import pandas as pd
import csv

df1 = pd.read_csv('/media/yuhua/yuhua_projects/enhProj/ENHData/filtered_mm10.csv',header = 0)

enh = df1[['chr','start','end','length']]
enh.index = list(range(1,len(enh)+1))

# enh_count0 = pd.read_table("/media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_GSR_counts.txt",index_col=0)
enh_count0 = pd.read_table("/media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_XW_counts.txt",index_col=0)
enh_count = enh_count0[enh_count0['gene_id'].str.contains("Enh")] 

# timelist_GSR=['2cell', '4cell', '8cell', 'EpiE65', 'ExeE65', 'ICM','MIIOocyte', 'morula', 'TE']
timelist_XW=['E2C','ICM','L2C','M4C','M8C','MIIOocyte','Zygote']

enh2 = pd.concat([enh,enh_count],axis=1)

def cal_EPM(df,timelist):
    df_co = df.copy()
    for i in timelist:
        df_co[i]=df[i]*1000/df['length']
    alstreads=[]
    for i in timelist:
        alstreads.append(sum(df_co[i]))
    k = 0
    df_co2=df_co.copy()
    for i in timelist:
        df_co2[i] = (df_co[i]/alstreads[k])*1000000
        k=k+1
    return df_co2

enh_EPM = cal_EPM(enh2,timelist_XW)
# enh_EPM.to_csv("/media/yuhua/yuhua_projects/enhProj/ENHData/enh_GSR_EPM.csv")
enh_EPM.to_csv("/media/yuhua/yuhua_projects/enhProj/ENHData/enh_XW_EPM.csv")
