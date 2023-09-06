## Map sequencing reads (BAM/fragments.tsv/.bed) to cPeaks

at map2cpeaks folder, just download it and use it.

### Version
Python (>=3.7)

### Requirments

```
anndata
scanpy
numpy
tqdm
joblib
```

### Method 1: Map the sequencing reads (fragments.tsv.gz) in each sample/cell to generate cell-by-cPeak matrix (.mtx/.h5ad)
 
```
usage:

1.


cd map2cpeaks


2. 

python main.py --fragment_path path/to/your_fragment.tsv.gz
               
--fragment_path, -f: the input file must be *.tsv.gz file, transfer fragment to cPeak reference
 
optional arguments:

 --help, -h:          show this help message
 --barcode_path, -b:  barcode file is given, the code will use the barcode in the file or the code will use all barcodes in the fragment.tsv.gz file
 --type_saved, -t:    output file is mtx file or h5ad former, Default to .mtx
 --output, -o:        output folder, Default to ./map2cpeaks_result
 --output_name:       name of output files, Default to cell-cpeaks.
 --num_cores, -n:     number of cores to use, Default to 10.
 --reference:         cPeak version, hg38 or hg19, Default to hg38.
```


### Method 2. Directly map the pre-identified features like peaks to cPeaks (NOT recommand)

**This is not a good idea.** It may lose information in the genomic regions which are not included in pre-identfied features. Also, for bulk ATAC-seq data, the quantification of each cPeak is inaccurate.

```
usage: python main.py [--bed_path feature.bed]

--bed_path, -bed: the input feature.bed file, for example, MACS2calledPeaks.bed.
```
