## Introduction

cPeaks serves as a unified reference for different cell types for scATAC-seq data, enhancing downstream analysis, particularly in cell annotation and rare cell type detection.

<img src=".\Tutorials\media\methods.png" alt="1" style="zoom:100%;" />

## cPeaks Download

We released cPeaks (.bed format) with [GRCh37/hg19](https://cloud.tsinghua.edu.cn/f/7b7664158dd7482c9a95/?dl=1) and [GRCh38/hg38](https://cloud.tsinghua.edu.cn/f/ff4591857f5d472d9401/?dl=1) version.

The basic information and properties of cPeaks can be found in [cPeaks_info.tsv](https://cloud.tsinghua.edu.cn/f/4422592d373948589dc4/).

## cPeaks Usage

### Quick Start Guide

cPeaks eliminates the need for the peak-calling step by providing a ready-to-use reference. If you are familiar with SnapATAC2, ArchR or Python, this section will help you quickly get started with cPeaks. However, if you have any questions, please refer to [detailed tutorials](#detail).

- SnapATAC2
  
  ```python
  # Read cPeaks file cpeaks_hg19.bed or cpeaks_hg38.bed
  cpeaks_path = 'YOUR_PATH/cpeaks_hg38.bed'
  with open(cpeaks_path) as cpeaks_file:
      cpeaks = cpeaks_file.read().strip().split('\n')
      cpeaks = [peak.split('\t')[0] + ':' + peak.split('\t')[1] + '-' + peak.split('\t')[2] for peak in cpeaks]
  # Set parameter use_rep to 'cpeaks' as mapping reference. 'data' is your AnnData object.
  data = snap.pp.make_peak_matrix(data, use_rep=cpeaks)
  ```

- ArchR
  
  ```r
  # Read cPeaks file cpeaks_hg19.bed or cpeaks_hg38.bed
  cpeaks <- read_table('YOUR_PATH/cpeaks_hg19.bed', col_names = F)
  cpeaks.gr <- GRanges(seqnames = cpeaks$X1, ranges = IRanges(cpeaks$X2, cpeaks$X3))
  # Set parameter 'features' to cpeaks.gr. 'proj' is your ArchRProject object.
  proj <- addFeatureMatrix(proj, features = cpeaks.gr, matrixName = 'FeatureMatrix')
  ```

- Run Python Script Manually
  
  ```bash
  git clone https://github.com/MengQiuchen/cPeaks.git
  cd cPeaks/map2cpeak
  python main.py --fragment_path PATH/to/YOUR_fragment.tsv.gz --output map2cpeaks_result
  ```
  
    After running `main.py` as above, you will get a matrix file `cell_cpeaks.mtx` and a text file `barcodes.txt` within a newly created folder named `map2cpeaks_result`.
    
    [Learn more about the arguments for main.py](https://mengqiuchen.github.io/cPeaks/Tutorials/#/?id=arguments).

### <a id="detail"></a>Comprehensive Guide

Please refer to [detailed tutorials](https://mengqiuchen.github.io/cPeaks/Tutorials/#/?id=_3-comprehensive-guide).

## Citation

Meng Q, Wu X, et al. Toward a generic feature set defined by consensus peaks as a consistent reference for ATAC-seq data. Preprint at *bioRxiv* [Preprint](https://doi.org/10.1101/2023.05.30.542889)(2023).

## Contact

Please reach out to Meng Qiuchen at mqc17@mails.tsinghua.edu.cn if you encounter any issues or have any recommendations.