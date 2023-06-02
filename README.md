# readme

We released cPeaks, GRCh37 (hg19) and GRCh38 (hg38) version.
cPeaks contained observed part and predicted part. If you want to download only observed or predicted part in cPeaks, you can find them in link: xxxxx


 Here is the BASH command? to map bulk ATAC-seq &scATAC-seq data to cpeaks reference

1. map the sequencing reads in each sample/cell to cPeak reference
# command: xxxx --xxx ---xx --input_type auto/fragment/bam

# parameters:
 --fragment_path(must-have): '-f', the input file is fragment.gz file, which is the output of fragment file in the 10x pipeline

 --barcode_path(optional): '-b', if you give the barcode file, the code will use the barcode in the file; if you don't give the barcode file, the code will use all barcode in the fragment.gz file

 --reference（optional)： hg38/hg19, default is hg38. unavailable

 --version：all(default), "observed" only use observed part as reference, "predicted" only use predicted part as reference. unavailable

 --output(default:./res): if you do note give the output(default: ./res)

 --type_saved(default:NULL): for bulk ATAC-seq, the type saved is .bed. for single-cell, the default output file is "mtx" file for single-cell, or you can assign the "type_saved" to "h5ad" former

 --input_type:auto:auto detection (default); fragments: xxxx, bam: unavailable bed: unavailable

 example command: python main.py -f test_fragment.tsv.gz -b test_bar.txt -o result


2. directly map the pre-identified features like peaks to cPeaks (not recommand)

 this is not a good idea. It may lose the open regions in the genomic regions which are not included in pre-identfied features.If the input is bulk ATAC-seq data, the quantification of each cPeaks is unbelievable.

# command: 
 --bed_path: '-bed', the input file is peak.bed file

 example command: python main.py -bed test.bed






