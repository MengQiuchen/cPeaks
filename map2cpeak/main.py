# python 
# 2023/08/14
# version 1.2
# @auther : Xinze Wu

import os
import sys
import gzip
from tqdm import tqdm
from scipy.io import mmread
from collections import OrderedDict
from joblib import Parallel, delayed
import numpy as np
import scanpy as sc
import anndata
import argparse
import pandas as pd


# map a bed to b bed, the output is a tuple of two numpy arrays: (indices, values)
def map_A_to_B(a_peaks):
    
    chromosomes = list(range(1, 23))
    chromosomes.extend(["X", "Y"])
    chrlist = ["chr" + str(chrom) for chrom in chromosomes]
    
    a_dict = {}
    for chrom in chrlist:
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
def frag2mtx(fragment_path,savepath,barcode_path,num_cores):
    '''
    fragment_path: path of fragment.gz file
    savepath: path to save mtx file
    barcode_path: path of barcode file, ,txt format,each line is a barcode; if None, use all barcode in fragment.gz file
    '''
    
    # get unique barcode

    barcodes = open(barcode_path, 'r').readlines()
    barcodes = [i.strip() for i in barcodes]
    print('num of barcode is {}'.format(len(barcodes)))


    # get the bed in each barcode
    print('start to get bed in each barcode')

    fragment_data = gzip.open(fragment_path, 'rb')
    
    bar2bed = OrderedDict()
    for barcode in barcodes:
        bar2bed[barcode] = []
    
    for line in tqdm(fragment_data):
        temp = line.decode().strip().split('\t')
        try:
            bar2bed[temp[3]].append((temp[0],int(temp[1]),int(temp[2])))
        except:
            continue
    print('num peak of random bar: ',len(bar2bed[barcodes[-1]]))
    print('num peak of random bar: ',len(bar2bed[barcodes[-2]]))
    print('num peak of random bar: ',len(bar2bed[barcodes[-3]]))


    print('mapping to cpeaks...')

    results = Parallel(n_jobs=num_cores, backend="multiprocessing")(
        delayed(map_A_to_B)(input_list) for input_list in tqdm(list(bar2bed.values()))
    )
    
    value_num = 0
    for res in results:
        value_num += len(res[0])
    
    print('the num of the value of mtx is: ',value_num)
    
    print('start to write mtx file')
    global mtx_path
    mtx_path= os.path.join(savepath,output_name+'.mtx')
    
    with open(mtx_path, 'w') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n')
        f.write(f'{len(b_peaks)} {len(barcodes)} {value_num}\n')
        
        for i,res in enumerate(results):
            for j in range(len(res[0])):
                f.write(f'{res[0][j]} {i+1} {res[1][j]}\n')
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

    chromosomes = list(range(1, 23))
    chromosomes.extend(["X", "Y"])
    chrlist = ["chr" + str(chrom) for chrom in chromosomes]
    
    a_dict = {}
    for chrom in chrlist:
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
    parser.add_argument("--num_cores",'-n', type=int, default= 10,help="num of cores to use")
    parser.add_argument("--bed_path", '-bed',type=str, default=None,help="bed file you called before")
    parser.add_argument("--reference",type=str, default='hg38',help="cpeaks version: hg38 or hg19")

    args = parser.parse_args()
   
    from time import time
    time_ = time()

    fragment_path = args.fragment_pathll
    barcode_path = args.barcode_path
    savepath = args.output
    save_type = args.type_saved
    num_cores = args.num_cores
    bed_path = args.bed_path
    reference = args.reference
    output_name = args.output_name

    if reference == 'hg38':
        cpeaks_path = 'cpeaks_hg38.bed'
    else:
        cpeaks_path = 'cpeaks_hg19.bed'

    # get cpeaks
    try:
        print('reading cpeaks')
        with open(cpeaks_path) as b_file:
            b_lines = b_file.readlines()
        b_peaks = []
        for line in tqdm(b_lines):
            
            fields = line.strip().split('\t')
            b_peaks.append((fields[0], int(fields[1]), int(fields[2])))  
    except:
        raise('there is sth wrong with cpeaks file')
    
    if not os.path.exists(savepath):
        os.makedirs(savepath)


    if bed_path is not None:
        print('attention: you use bed file you called before')
        map_bed_to_bed(bed_path, b_peaks,os.path.join(savepath,'res.bed'))
        # 退出
        sys.exit(0)

    # fragment_path not end with 'tsv.gz' and not None
    if fragment_path[-6:] != 'tsv.gz':
        print(fragment_path[-6:])
        raise('fragment_path not end with tsv.gz')
    
    if barcode_path is None:
        
        barcode_path = os.path.join(savepath,'barcodes.txt')
       
        print('start to get barcode')
        fragment_data = gzip.open(fragment_path, 'rb')

        barcodes = set()
        for line in tqdm(fragment_data):
            line = line.decode()
            barcodes.add(line.split('\t')[3])
        print('save {} barcodes to {}'.format(len(barcodes),barcode_path))
        # write barcode to file
        
        with open(os.path.join(savepath,'barcodes.txt'), 'w') as f:
            for barcode in barcodes:
                f.write(barcode+'\n')


    if save_type == 'mtx':
        frag2mtx(fragment_path,savepath,barcode_path,num_cores)
    else:
        
        frag2mtx(fragment_path,savepath,barcode_path,num_cores)
        mtx2h5ad(os.path.join(savepath,output_name+'.mtx'), cpeaks_path, savepath)

    print('use time: ',time()-time_)
    
