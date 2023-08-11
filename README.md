## Introduction

Chromatin accessibility profiling methods such as assay for transposase-accessible chromatin using sequencing (ATAC-seq) have been promoting the identification of gene regulatory elements and the characterization of epigenetic landscapes. Unlike gene expression data, there is no consistent reference for chromatin accessibility data, which hinders large-scale integration analysis. By analyzing many more than 1000 ATAC-seq samples and 100 scATAC-seq samples, we found that cells share the same set of potential open regions. We thus proposed a reference called consensus peaks (cPeaks) to represent open regions across different cell types, and developed a deep-learning model to predict all potential open regions in the human genome. We showed that cPeaks can be regarded as a new set of epigenomic elements in the human genome, and using cPeaks can increase the performance of cell annotations and facilitate the discovery of rare cell types. cPeaks also performed well in analyzing dynamic biological processes and diseases. cPeaks can serve as a general reference for epigenetic studies, much like the reference genome for genomic studies, making the research faster, more accurate, and more scalable.

## Citation

Meng Q, Wu X, Li C, et al. The full set of potential open regions (PORs) in the human genome defined by consensus peaks of ATAC-seq data. [doi:10.1101/2023.05.30.542889](https://doi.org/10.1101/2023.05.30.542889)

## cPeaks

We released cPeaks (.bed format) with [GRCh37/hg19](https://github.com/MengQiuchen/cPeaks/blob/main/cPeaks_hg19.bed) and [GRCh38/hg38](https://github.com/MengQiuchen/cPeaks/blob/main/cPeaks_hg38.bed) version.

The basic information and properties of cPeaks can be found in [cPeaks_info.tsv](https://cloud.tsinghua.edu.cn/f/8fa73ab35ec049cf8358/).
- ID: The ID of this cPeak. We employed a 14-character string as the unique ID for each cPeak. cPeak ID begins with POR for potential open regions, then HS for humans, and ends with a nine-digital number. For example, the first cPeak was encoded as PORHS000000001.
- chr: The chromosome of this cPeak.
- source: The source of this cPeak, “observed” or “predicted”.
- start & end: start and end positions of this cPeak.
- housekeeping: the housekeeping status, "true" or "false".
- shapePattern: the shape pattern, "well-positioned", "one-side", "weakly-positioned" or "else".
- inferredElements: the inferred regulatory elements, "TSS", "Promoter", "Enhancer", "CTCF" or "TES". A cPeak may act as multiple inferred regulatory elements, which were split by ",". 

## Map sequencing reads (BAM/fragments.tsv/.bed) to cPeaks

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

### Method 1: Map the sequencing reads (BAM/fragments.tsv) in each sample/cell to generate cell-by-cPeak matrix (.mtx/.h5ad)
 
```
usage: python main.py [--fragment_path your_fragment.tsv.gz]
                      [--barcode_path your_bar.txt]
                      [--output your_output_folder]

 --fragment_path, -f: the input file is fragment.gz file, transfer fragment to cPeak reference
 
optional arguments:

 --help,-h:           show this help message
 --barcode_path, -b:  barcode file given, the code will use the barcode in the file or the code will use all barcodes in the fragment.gz file
 --type_saved，-t:    output file is mtx file or h5ad former, Default to .mtx
 --output, -o:        output folder, Default to ./res
 --reference:         cPeak version, hg38 or hg19, Default to hg38.

```
note: The observed part only use observed data as reference, the predicted part only use predicted cPeaks by DeepCNN as reference, and all is the combination of both.

### Method 2. Directly map the pre-identified features like peaks to cPeaks (NOT recommand)

**This is not a good idea.** It may lose information in the genomic regions which are not included in pre-identfied features. Also, for bulk ATAC-seq data, the quantification of each cPeak is inaccurate.

```
usage: python main.py [--bed_path feature.bed]

--bed_path, -bed: the input feature.bed file, for example: MACS2calledPeaks.bed.
```

## Contact us
mqc17@mails.tsinghua.edu.cn
