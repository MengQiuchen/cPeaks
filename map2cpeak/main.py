# python 
# 2024-3-8
# @auther : Xinze Wu

import os
import sys
import gzip
# import anndata
import argparse
# import pandas as pd
import numpy as np
# import scanpy as sc
from tqdm import tqdm
# from scipy.io import mmread
from collections import Counter,OrderedDict


# def function to trans fragment file to cpeaks referenced mtx
def frag2mtx(fragment_path,savepath,barcode_path):
    '''
    fragment_path: path of fragment.gz file
    savepath: path to save mtx file
    barcode_path: path of barcode file, ,txt format,each line is a barcode; if None, use all barcode in fragment.gz file
    '''
    bar2idx = OrderedDict()
    
    if barcode_path == None:


        # 打开压缩文件并且逐行读取
        with gzip.open(fragment_path, 'rt') as file:
            #          chr     start   end          barcode
            # array([['chr1', 181500, 181531, 'AAACAGCCAACCCTAA-1'],
            #        ['chr1', 629913, 629992, 'AAACAGCCAACCCTAA-1'],
            #        ['chr1', 629914, 629985, 'AAACAGCCAACCCTAA-1'],
            #        ['chr1', 629914, 629986, 'AAACAGCCAACCCTAA-1'],
            #        ['chr1', 629939, 629986, 'AAACAGCCAACCCTAA-1'],
            #        ['chr1', 633935, 634098, 'AAACAGCCAACCCTAA-1'],
            #        ['chr1', 633957, 634087, 'AAACAGCCAACCCTAA-1'],
            #        ... ])

            for line in tqdm(file):
                tmp = line.strip().split('\t')
                if len(tmp) < 5 or tmp[0] not in chr_list:
                    continue
                if tmp[3] not in bar2idx:
                    bar2idx[tmp[3]] = []

                idx_set = set()
                for i in range(int(tmp[1]), int(tmp[2]), 59):
                    try:
                        if dic_chr[tmp[0]][i]!=-1:
                            idx_set.add(dic_chr[tmp[0]][i])
                    except:
                        continue
                try:
                    if dic_chr[tmp[0]][int(tmp[2])-1]!=-1:
                        idx_set.add(dic_chr[tmp[0]][int(tmp[2])-1])
                except:
                    continue
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
                    try:
                        if dic_chr[tmp[0]][i]!=-1:
                            idx_set.add(dic_chr[tmp[0]][i])
                    except:
                        continue
                try:
                    if dic_chr[tmp[0]][int(tmp[2])-1]!=-1:
                        idx_set.add(dic_chr[tmp[0]][int(tmp[2])-1])
                except:
                    continue
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

# def mtx2h5ad(mtx_path, barcode_path,savepath):

#     '''
#     mtx_path: path of mtx file
#     savepath: path to save h5ad file
#     '''
#     print('read mtx')
#     mtx = mmread(mtx_path).tocsr()
    
#     cpeaks_index = [f'{i[0]}:{i[1]}-{i[2]}' for i in b_peaks]
    
#     print('read barcodes')
#     barcodes = [i.strip() for i in open(barcode_path).readlines()]
    
#     adata = sc.AnnData(X = mtx.T, obs =  {'barcodes':barcodes},var = {'cpeaks':cpeaks_index},dtype=mtx.dtype)
#     adata.obs_names = barcodes
    
#     adata.write(os.path.join(savepath,output_name+'.h5ad'))
    
#     print('write to h5ad done!')

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

def cpeak_init(mode,reference):

    
    # get cpeaks
    # print('load cpeaks...')

    # hg38
    # chr1    248956422
    # chr2    242193529
    # chr3    198295559
    # chr4    190214555
    # chr5    181538259
    # chr6    170805979
    # chr7    159345973
    # chr8    145138636
    # chr9    138394717
    # chr10   133797422
    # chr11   135086622
    # chr12   133275309
    # chr13   114364328
    # chr14   107043718
    # chr15   101991189
    # chr16   90338345
    # chr17   83257441
    # chr18   80373285
    # chr19   58617616
    # chr20   64444167
    # chr21   46709983
    # chr22   50818468
    # chrX    156040895
    # chrY    57227415

    # hg19
    # chr1    249250621
    # chr2    243199373
    # chr3    198022430
    # chr4    191154276
    # chr5    180915260
    # chr6    171115067
    # chr7    159138663
    # chr8    146364022
    # chr9    141213431
    # chr10   135534747
    # chr11   135006516
    # chr12   133851895
    # chr13   115169878
    # chr14   107349540
    # chr15   102531392
    # chr16   90354753
    # chr17   81195210
    # chr18   78077248
    # chr19   59128983
    # chr20   63025520
    # chr21   48129895
    # chr22   51304566
    # chrX    155270560
    # chrY    59373566


    dic_chr = {}
    if mode == 'performance':
        if reference == 'hg38':
            dic_chr['chr1'] = np.full(248956422, -1,dtype=int)
            dic_chr['chr2'] = np.full(242193529, -1,dtype=int)
            dic_chr['chr3'] = np.full(198295559, -1,dtype=int)
            dic_chr['chr4'] = np.full(190214555, -1,dtype=int)
            dic_chr['chr5'] = np.full(181538259, -1,dtype=int)
            dic_chr['chr6'] = np.full(170805979, -1,dtype=int)
            dic_chr['chr7'] = np.full(159345973, -1,dtype=int)
            dic_chr['chr8'] = np.full(145138636, -1, dtype=int)
            dic_chr['chr9'] = np.full(138394717, -1, dtype=int)
            dic_chr['chr10'] = np.full(133797422, -1, dtype=int)
            dic_chr['chr11'] = np.full(135086622, -1, dtype=int)
            dic_chr['chr12'] = np.full(133275309, -1, dtype=int)
            dic_chr['chr13'] = np.full(114364328, -1, dtype=int)
            dic_chr['chr14'] = np.full(107043718, -1, dtype=int)
            dic_chr['chr15'] = np.full(101991189, -1, dtype=int)
            dic_chr['chr16'] = np.full(90338345, -1, dtype=int)
            dic_chr['chr17'] = np.full(83257441, -1, dtype=int)
            dic_chr['chr18'] = np.full(80373285, -1, dtype=int)
            dic_chr['chr19'] = np.full(58617616, -1, dtype=int)
            dic_chr['chr20'] = np.full(64444167, -1, dtype=int)
            dic_chr['chr21'] = np.full(46709983, -1, dtype=int)
            dic_chr['chr22'] = np.full(50818468, -1, dtype=int)
            dic_chr['chrX'] = np.full(156040895, -1, dtype=int)
            dic_chr['chrY'] = np.full(57227415, -1, dtype=int)
        elif reference == 'hg19':
            dic_chr['chr1'] = np.full(249250621, -1, dtype=int)
            dic_chr['chr2'] = np.full(243199373, -1, dtype=int)
            dic_chr['chr3'] = np.full(198022430, -1, dtype=int)
            dic_chr['chr4'] = np.full(191154276, -1, dtype=int)
            dic_chr['chr5'] = np.full(180915260, -1, dtype=int)
            dic_chr['chr6'] = np.full(171115067, -1, dtype=int)
            dic_chr['chr7'] = np.full(159138663, -1, dtype=int)
            dic_chr['chr8'] = np.full(146364022, -1, dtype=int)
            dic_chr['chr9'] = np.full(141213431, -1, dtype=int)
            dic_chr['chr10'] = np.full(135534747, -1, dtype=int)
            dic_chr['chr11'] = np.full(135006516, -1, dtype=int)
            dic_chr['chr12'] = np.full(133851895, -1, dtype=int)
            dic_chr['chr13'] = np.full(115169878, -1, dtype=int)
            dic_chr['chr14'] = np.full(107349540, -1, dtype=int)
            dic_chr['chr15'] = np.full(102531392, -1, dtype=int)
            dic_chr['chr16'] = np.full(90354753, -1, dtype=int)
            dic_chr['chr17'] = np.full(81195210, -1, dtype=int)
            dic_chr['chr18'] = np.full(78077248, -1, dtype=int)
            dic_chr['chr19'] = np.full(59128983, -1, dtype=int)
            dic_chr['chr20'] = np.full(63025520, -1, dtype=int)
            dic_chr['chr21'] = np.full(48129895, -1, dtype=int)
            dic_chr['chr22'] = np.full(51304566, -1, dtype=int)
            dic_chr['chrX'] = np.full(155270560, -1, dtype=int)
            dic_chr['chrY'] = np.full(59373566, -1, dtype=int)

        else:
            raise ValueError('reference must be hg38 or hg19')
    else:
        int32 = np.int32
        if reference == 'hg38':
            dic_chr['chr1'] = np.full(248956422, -1,dtype=int32)
            dic_chr['chr2'] = np.full(242193529, -1,dtype=int32)
            dic_chr['chr3'] = np.full(198295559, -1,dtype=int32)
            dic_chr['chr4'] = np.full(190214555, -1,dtype=int32)
            dic_chr['chr5'] = np.full(181538259, -1,dtype=int32)
            dic_chr['chr6'] = np.full(170805979, -1,dtype=int32)
            dic_chr['chr7'] = np.full(159345973, -1,dtype=int32)
            dic_chr['chr8'] = np.full(145138636, -1,dtype=int32)
            dic_chr['chr9'] = np.full(138394717, -1,dtype=int32)
            dic_chr['chr10'] = np.full(133797422, -1,dtype=int32)
            dic_chr['chr11'] = np.full(135086622, -1,dtype=int32)
            dic_chr['chr12'] = np.full(133275309, -1,dtype=int32)
            dic_chr['chr13'] = np.full(114364328, -1,dtype=int32)
            dic_chr['chr14'] = np.full(107043718, -1,dtype=int32)
            dic_chr['chr15'] = np.full(101991189, -1,dtype=int32)
            dic_chr['chr16'] = np.full(90338345, -1,dtype=int32)
            dic_chr['chr17'] = np.full(83257441, -1,dtype=int32)
            dic_chr['chr18'] = np.full(80373285, -1,dtype=int32)
            dic_chr['chr19'] = np.full(58617616, -1,dtype=int32)
            dic_chr['chr20'] = np.full(64444167, -1,dtype=int32)
            dic_chr['chr21'] = np.full(46709983, -1,dtype=int32)
            dic_chr['chr22'] = np.full(50818468, -1,dtype=int32)
            dic_chr['chrX'] = np.full(156040895, -1,dtype=int32)
            dic_chr['chrY'] = np.full(57227415, -1,dtype=int32)

        elif reference == 'hg19':
            dic_chr['chr1'] = np.full(249250621, -1,dtype=int32)
            dic_chr['chr2'] = np.full(243199373, -1,dtype=int32)
            dic_chr['chr3'] = np.full(198022430, -1,dtype=int32)
            dic_chr['chr4'] = np.full(191154276, -1,dtype=int32)
            dic_chr['chr5'] = np.full(180915260, -1,dtype=int32)
            dic_chr['chr6'] = np.full(171115067, -1,dtype=int32)
            dic_chr['chr7'] = np.full(159138663, -1,dtype=int32)
            dic_chr['chr8'] = np.full(146364022, -1,dtype=int32)
            dic_chr['chr9'] = np.full(141213431, -1,dtype=int32)
            dic_chr['chr10'] = np.full(135534747, -1,dtype=int32)
            dic_chr['chr11'] = np.full(135006516, -1,dtype=int32)
            dic_chr['chr12'] = np.full(133851895, -1,dtype=int32)
            dic_chr['chr13'] = np.full(115169878, -1,dtype=int32)
            dic_chr['chr14'] = np.full(107349540, -1,dtype=int32)
            dic_chr['chr15'] = np.full(102531392, -1,dtype=int32)
            dic_chr['chr16'] = np.full(90354753, -1,dtype=int32)
            dic_chr['chr17'] = np.full(81195210, -1,dtype=int32)
            dic_chr['chr18'] = np.full(78077248, -1,dtype=int32)
            dic_chr['chr19'] = np.full(59128983, -1,dtype=int32)
            dic_chr['chr20'] = np.full(63025520, -1,dtype=int32)
            dic_chr['chr21'] = np.full(48129895, -1,dtype=int32)
            dic_chr['chr22'] = np.full(51304566, -1,dtype=int32)
            dic_chr['chrX'] = np.full(155270560, -1,dtype=int32)
            dic_chr['chrY'] = np.full(59373566, -1,dtype=int32)
        else:
            raise ValueError('reference must be hg38 or hg19')
    return dic_chr
 



    
  

# main 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fragment_path", '-f',type=str, default=None,help="path of mtx file")
    parser.add_argument("--barcode_path", '-b',type=str, default=None,help="path of barcode file")
    parser.add_argument("--output", '-o',type=str, default='map2cpeaks_result',help="path to save files")
    parser.add_argument("--output_name",type=str, default='cell-cpeaks',help="name of output files e.g. cell-cpeaks")
    parser.add_argument("--mode",type=str, default='performance',help="performance or normal, performance mode enabled: More memory will be used")
    parser.add_argument("--bed_path", '-bed',type=str, default=None,help="bed file you called before")
    parser.add_argument("--reference",type=str, default='hg38',help="cpeaks version: hg38 or hg19")

    args = parser.parse_args()
    
   
    from time import time
    time_ = time()

    fragment_path = args.fragment_path
    barcode_path = args.barcode_path
    savepath = args.output
    mode = args.mode
    bed_path = args.bed_path
    reference = args.reference
    output_name = args.output_name
    
    # fragment_path not end with 'tsv.gz' and not None
    try:
        with gzip.open(fragment_path, 'rt') as file:
            for line in file:
                if line[0]=='#' or line[0]=='\n':
                    continue
                tmp = line.strip().split('\t')
                if len(tmp)<3:
                    raise ValueError('tsv file is empty or the separator is not \\t')
                break
    except:
        raise ValueError('Input file is not a valid tsv.gz file')
    
    if barcode_path is None:
        
        print("using all barcodes in the fragment file")
    try:
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        print('file will save to ',savepath)
    except:
        raise ValueError('savepath is not a valid path or you do not have the permission to create the folder')
            
    chr_list = ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']

    if reference == 'hg38':
        cpeaks_path = 'cpeaks_hg38.bed.gz'
    elif reference == 'hg19':
        cpeaks_path = 'cpeaks_hg19.bed.gz'
    else:
        raise ValueError('reference must be hg38 or hg19')

    print('using reference:',cpeaks_path)

    dic_chr = cpeak_init(mode,reference)
   
    try:
        b_peaks = []
        print("read and load cpeaks...")
        with gzip.open(cpeaks_path,'rt') as b_file:
            b_lines = b_file.readlines()
            
        for idx,line in tqdm(enumerate(b_lines)):
            fields = line.strip().split('\t')
            b_peaks.append((fields[0], int(fields[1]), int(fields[2])))
            dic_chr[fields[0]][int(fields[1]):(int(fields[2])+1)] = idx
            
        print('load cpeak finished!')
             
    except:
        raise ValueError('there is sth wrong with cpeaks file, plz send email to wxz22@mails.tsinghua.edu.cn')
    
    
    if bed_path is not None:
        print('attention: you use bed file you called before')
        map_bed_to_bed(bed_path, b_peaks,os.path.join(savepath,'res.bed'))
        print("finished!")
        sys.exit(0)
    
    # if save_type == 'mtx':
    frag2mtx(fragment_path,savepath,barcode_path)
    # elif save_type == 'h5ad':
        
    #     frag2mtx(fragment_path,savepath,barcode_path)
    #     if barcode_path is None:
    #         mtx2h5ad(os.path.join(savepath,output_name+'.mtx'),os.path.join(savepath,'barcodes.txt'), savepath)
    #     else:
    #         mtx2h5ad(os.path.join(savepath,output_name+'.mtx'),barcode_path, savepath)
            
            
    # else:
    #     raise('--type_saved (-t) must be mtx or h5ad')
        

    print('use time: ',round(time()-time_),"s")
    
