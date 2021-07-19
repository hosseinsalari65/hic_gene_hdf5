import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from time import time
from matplotlib import cm
chromosome=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr1\
5','chr16','chr17','chr18','chr19','chrX']

c=0
ch = chromosome[c]
#lg0 = 2000
#df_genes = pd.read_csv(ch+"/gene_list_RNAseq_"+ch+"_h"+str(int(lg0/1000))+"kb.txt",sep="\t")

#start0 =df_genes.start[0]
#end0   = df_genes.end[0]
#print(df_genes.length[0],df_genes.strand[0],df_genes.RNAseq_fpkm[0])

#start0 = 0
#end0 = 10000000

res = 10000
val = 'resolutions'
normal = 'KR'
#start_bin = int(start0/res)
#end_bin = int(end0/res)
#length_matrix = end_bin-start_bin+1
#print(start_bin,end_bin,length_matrix)

#HiC2 = np.zeros((length_matrix,length_matrix))
#HiC3 = np.zeros((length_matrix,length_matrix))
#HiC4 = np.zeros((length_matrix,length_matrix))
with h5py.File("../mouse_ES_micro-c/4DNFINNZDDXV.mcool", "r") as h5f:
    # Get a h5py dataset object
    t0 = time()
    print(list(h5f[val]))
    chromid = h5f[val][str(res)]['bins']['chrom']
    bin_start = h5f[val][str(res)]['bins']['start']
    bin_end = h5f[val][str(res)]['bins']['end']
    bin_weigth1 = h5f[val][str(res)]['bins']['KR']
    bin_weigth2 = h5f[val][str(res)]['bins']['VC']
    bin_weigth3 = h5f[val][str(res)]['bins']['VC_SQRT']
    chroms_name = h5f[val][str(res)]['chroms']['name']
    chroms_length = h5f[val][str(res)]['chroms']['length']
    chrom_offset =h5f[val][str(res)]['indexes']['chrom_offset']
    bin1_offset = h5f[val][str(res)]['indexes']['bin1_offset']
    for i in range(len(chroms_name)):
        if (chroms_name[i].decode('UTF-8')==ch):
            chID = i
            break
    start_bin = chrom_offset[chID]
    end_bin = chrom_offset[chID+1]
    print(start_bin,end_bin,end_bin-start_bin)
    arr_start = bin1_offset[start_bin]
    arr_end = bin1_offset[end_bin]+1
    length_chr = end_bin - start_bin +1
    HiC1 = np.zeros((length_chr,length_chr))
    print(arr_start,arr_end,length_chr)
    bin1_t = h5f[val][str(res)]['pixels']['bin1_id']
    bin2_t = h5f[val][str(res)]['pixels']['bin2_id']
    count_t=h5f[val][str(res)]['pixels']['count']
    df = pd.DataFrame({'bin1':bin1_t[arr_start:arr_end],'bin2':bin2_t[arr_start:arr_end],'raw':count_t[arr_start:arr_end]})
    print(len(df))
    sub_df = df[(df['bin2']<=end_bin)]
#    sub_df['wt1'] = sub_df['bin1_t']
    sub_df['wt1'] = sub_df['bin1'].apply(lambda x: bin_weigth1[x])
#    sub_df['wt2'] = sub_df['bin2_t']
    sub_df['wt2'] = sub_df['bin1'].apply(lambda x: bin_weigth1[x])
    sub_df['norm'] = sub_df['wt1']*sub_df['wt2']*sub_df['raw']
    sub_df = sub_df.drop(['wt1','wt2'], 1)
    print(sub_df)
    print(sub_df.shape)
#    sub_df = sub_df[(sub_df['bin2_t']>arr_start)&(sub_df['bin2_t']<arr_end)]
#    print(len(sub_df))
    l = len(sub_df)
#    for i in range(arr_start,arr_end+1):
    for i in range(l):
        idx = sub_df.index[i]
        bin1 = sub_df.bin1_t[idx]
        bin2 = sub_df.bin2_t[idx]
        count= sub_df.count_t[idx]
        wt1 = bin_weigth1[bin1]*bin_weigth1[bin2]*count
        HiC1[bin1-start_bin,bin2-start_bin] = (wt1)
        HiC1[bin2-start_bin,bin1-start_bin] = (wt1)
#        wt2 = bin_weigth2[bin1]*bin_weigth2[bin2]*count
#        HiC2[bin1-start_bin,bin2-start_bin] = (wt2)
#        HiC2[bin2-start_bin,bin1-start_bin] = (wt2)
#        wt3 = bin_weigth3[bin1]*bin_weigth3[bin2]*count
#        HiC3[bin1-start_bin,bin2-start_bin] = (wt3)
#        HiC3[bin2-start_bin,bin1-start_bin] = (wt3)
#        HiC4[bin1-start_bin,bin2-start_bin] = (count)
#        HiC4[bin2-start_bin,bin1-start_bin] = (count)
            #print(bin2-bin1,wt)
            #contact[bin2-bin1].append()
        #print(bin1,bin2,count)
    
    
    print(time()-t0)
    #plt.subplot(2,2,1)
    #plt.imshow((HiC1), vmin=-5, vmax=5,cmap='YlOrRd')
    #plt.imshow(HiC1,cmap='YlOrRd')
    plt.imshow(np.log2(HiC1),cmap='YlOrRd')
    plt.title(normal+', resolution='+str(res/1000)+'kb')
    plt.xlabel(ch+':'+str(start0/1000000)+'-'+str(end0/1000000)+'Mb')
    plt.ylabel(ch+':'+str(start0/1000000)+'-'+str(end0/1000000)+'Mb')
    plt.xticks([])
    plt.yticks([])
    #plt.subplot(2,2,2)
    #plt.imshow(np.log2(HiC2),cmap='YlOrRd')
    #plt.title('VC')
    #plt.subplot(2,2,3)
    #plt.imshow(np.log2(HiC3),cmap='YlOrRd')
    #plt.title('VC_SQRT')
    #plt.subplot(2,2,4)
    #plt.imshow(np.log2(HiC4),cmap='YlOrRd')
    #plt.title('raw')
    cbar = plt.colorbar()
    cbar.set_label('contact frequency (Log2)')
    plt.show()
