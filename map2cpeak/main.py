# python 
# 2023-09-05
# version 2.0
# @auther : Xinze Wus

import os
import sys
import gzip
import anndata
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
from tqdm import tqdm
from scipy.io import mmread
from collections import Counter,OrderedDict



# map a bed to b bed, the output is a tuple of two numpy arrays: (indices, values)
def map_A_to_B(a_peaks):
    
    a_dict = {}
    for chrom in chr_list:
        a_dict[chrom] = []
    
    for peak in a_peaks:
        chrom = peak[0]
        start = peak[1]
        end = peak[2]
        if chrom in a_dict:
            a_dict[chrom].append((start, end))
    
    overlap_counts = np.zeros(len(b_peaks), dtype=int)
    for i, peak in enumerate(b_peaks):
        chrom = peak[0]
        start = peak[1]
        end = peak[2]
        if chrom not in a_dict:
            continue
        
        for a_peak in a_dict[chrom]:
            if a_peak[0] < end and start < a_peak[1]:
                overlap_counts[i] += 1
          
    non_zero_indices = np.nonzero(overlap_counts)[0]
    non_zero_values = overlap_counts[non_zero_indices]
    
    return (non_zero_indices+1, non_zero_values)

# def function to trans fragment file to cpeaks referenced mtx
def frag2mtx(fragment_path,savepath,barcode_path):
    '''
    fragment_path: path of fragment.gz file
    savepath: path to save mtx file
    barcode_path: path of barcode file, ,txt format,each line is a barcode; if None, use all barcode in fragment.gz file
    '''
    bar2idx = OrderedDict()
    
    if barcode_path == None:
        
        with gzip.open(fragment_path, 'rt') as file:
            for line in tqdm(file):
                tmp = line.strip().split('\t')
                if len(tmp) < 5 or tmp[0] not in chr_list:
                    continue
                if tmp[3] not in bar2idx:
                    bar2idx[tmp[3]] = []

                idx_set = set()
                for i in range(int(tmp[1]), int(tmp[2]), 59):
                    if i in dic_chr[tmp[0]].keys():
                        idx_set.add(dic_chr[tmp[0]][i])

                if int(tmp[2])-1 in dic_chr[tmp[0]]:
                    idx_set.add(dic_chr[tmp[0]][int(tmp[2])-1])
                bar2idx[tmp[3]] += list(idx_set)

        value_num = 0
        bar2dic = OrderedDict()
        for i,j in bar2idx.items():
            counter = Counter(j)
            bar2dic[i] = counter
            value_num+=len(counter)

        fb = open(os.path.join(savepath,'barcodes.txt'),'w')
        fm = open(os.path.join(savepath,output_name+'.mtx'),'w')


        fm.write('%%MatrixMarket matrix coordinate real general\n')
        fm.write(f'{len(b_peaks)} {len(bar2idx)} {value_num}\n')

        for idx,(i,dic) in enumerate(bar2dic.items()):
            fb.write(i+'\n')
            for c,v in dic.items():
                fm.write(f'{c+1} {idx+1} {v}\n')

        fb.close()
        fm.close()
    else:
        barcodes = [i.strip() for i in open(barcode_path).readlines()]
        for barcode in barcodes:
            bar2idx[barcode] = []
            
        with gzip.open(fragment_path, 'rt') as file:
            for line in tqdm(file):
                tmp = line.strip().split('\t')
                if len(tmp) < 5 or tmp[0] not in chr_list or tmp[3] not in bar2idx:
                    continue

                idx_set = set()
                for i in range(int(tmp[1]), int(tmp[2]), 59):
                    if i in dic_chr[tmp[0]].keys():
                        idx_set.add(dic_chr[tmp[0]][i])

                if int(tmp[2])-1 in dic_chr[tmp[0]]:
                    idx_set.add(dic_chr[tmp[0]][int(tmp[2])-1])
                bar2idx[tmp[3]] += list(idx_set)
        value_num = 0
        bar2dic = OrderedDict()
        for i,j in bar2idx.items():
            counter = Counter(j)
            bar2dic[i] = counter
            value_num+=len(counter)

        fm = open(os.path.join(savepath,output_name+'.mtx'),'w')

        fm.write('%%MatrixMarket matrix coordinate real general\n')
        fm.write(f'{len(b_peaks)} {len(bar2idx)} {value_num}\n')

        for idx,(i,dic) in enumerate(bar2dic.items()):
            for c,v in dic.items():
                fm.write(f'{c+1} {idx+1} {v}\n')

        fm.close()
            
    print('fragment to mtx done!')


# mtx, cpeaks, barcode to h5ad
def mtx2h5ad(mtx_path, cpeaks_path, savepath):
    '''
    mtx_path: path of mtx file
    cpeaks_path: path of cpeaks file
    savepath: path to save h5ad file
    '''
    print('read mtx')
    mtx = mmread(mtx_path).tocsr()
    
    cpeaks_index = [f'{i[0]}:{i[1]}-{i[2]}' for i in b_peaks]
    
    print('read barcodes')
    barcodes = [i.strip() for i in open(barcode_path).readlines()]
    
    adata = sc.AnnData(X = mtx.T, obs =  {'barcodes':barcodes},var = {'cpeaks':cpeaks_index},dtype=mtx.dtype)
    adata.obs_names = barcodes
    
    adata.write(os.path.join(savepath,output_name+'.h5ad'))
    
    print('write to h5ad done!')

# map bed to bed 
def map_bed_to_bed(a_bed_path, b_peaks,savepath):

    """
    a_bed_path: path of a bed file
    b_bed: path of cpeak bed
    savepath: path to save result
    """ 

    a_peaks = open(a_bed_path, 'r').readlines()
    a_peaks = [i.strip().split('\t') for i in a_peaks]

    
    a_dict = {}
    for chrom in chr_list:
        a_dict[chrom] = []
    
    for peak in a_peaks:
        chrom = peak[0]
        start = int(peak[1])
        end = int(peak[2])
        if chrom in a_dict:
            a_dict[chrom].append((start, end))
    
    overlap_counts = np.zeros(len(b_peaks), dtype=int)
    for i, peak in tqdm(enumerate(b_peaks)):
        chrom = peak[0]
        start = int(peak[1])
        end = int(peak[2])
        if chrom not in a_dict:
            continue
        
        for a_peak in a_dict[chrom]:
            if a_peak[0] < end and start < a_peak[1]:
                overlap_counts[i] += 1
    # get the index of non-zero value
    non_zero_indices = np.nonzero(overlap_counts)[0]
    # map to cpeaks and write to bed
    with open(savepath, 'w') as f:
        for i in non_zero_indices:
            f.write(b_peaks[i][0]+'\t'+str(b_peaks[i][1])+'\t'+str(b_peaks[i][2])+'\n')
    print('done') 

# main 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fragment_path", '-f',type=str, default=None,help="path of mtx file")
    parser.add_argument("--barcode_path", '-b',type=str, default=None,help="path of barcode file")
    parser.add_argument("--output", '-o',type=str, default='./map2cpeaks_result',help="path to save files")
    parser.add_argument("--output_name",type=str, default='cell-cpeaks',help="name of output files e.g. cell-cpeaks")
    parser.add_argument("--type_saved", '-t',type=str, default='mtx',help="save type, h5ad or mtx")
    parser.add_argument("--bed_path", '-bed',type=str, default=None,help="bed file you called before")
    parser.add_argument("--reference",type=str, default='hg38',help="cpeaks version: hg38 or hg19")

    args = parser.parse_args()
    
   
    from time import time
    time_ = time()

    fragment_path = args.fragment_path
    barcode_path = args.barcode_path
    savepath = args.output
    save_type = args.type_saved
    bed_path = args.bed_path
    reference = args.reference
    output_name = args.output_name
    
    # fragment_path not end with 'tsv.gz' and not None
    if fragment_path == None or fragment_path[-6:] != 'tsv.gz':
        raise ValueError('fragment_path does not end with tsv.gz')
    
    if barcode_path is None:
        
        print("using all barcodes in the fragment file")
        
    
    chr_list = ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']


    if reference == 'hg38':
        cpeaks_path = 'cpeaks_hg38.bed.gz'
    else:
        cpeaks_path = 'cpeaks_hg19.bed.gz'

    # get cpeaks
    dic_chr = OrderedDict()
    b_peaks = []
    for i in range(1,23):
        dic_chr[f'chr{i}'] = OrderedDict()
    dic_chr['chrX'] = OrderedDict()
    dic_chr['chrY'] = OrderedDict()
    
    try:
        print("read and load cpeaks... 3 min is need...",end = ' ')
        with gzip.open(cpeaks_path,'rt') as b_file:
            b_lines = b_file.readlines()
            
        for idx,line in enumerate(b_lines):
            fields = line.strip().split('\t')
            b_peaks.append((fields[0], int(fields[1]), int(fields[2]))) 
            for j in range(int(fields[1]),int(fields[2])):
                dic_chr[fields[0]][j] = idx
        print('done')
             
    except:
        raise('there is sth wrong with cpeaks file, plz send email to wxz22@mails.tsinghua.edu.cn')
    
    if not os.path.exists(savepath):
        os.makedirs(savepath)

    if bed_path is not None:
        print('attention: you use bed file you called before')
        map_bed_to_bed(bed_path, b_peaks,os.path.join(savepath,'res.bed'))
        print("finished!")
        sys.exit(0)

    
    if save_type == 'mtx':
        frag2mtx(fragment_path,savepath,barcode_path)
    else:
        
        frag2mtx(fragment_path,savepath,barcode_path,num_cores)
        mtx2h5ad(os.path.join(savepath,output_name+'.mtx'), cpeaks_path, savepath)

    print('use time: ',round(time()-time_),"s")
    
