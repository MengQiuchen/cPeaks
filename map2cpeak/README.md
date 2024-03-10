# Map Sequencing Reads to cPeaks Using Map2cpeak

## Getting Started
Navigate to the `map2cpeaks` directory for the setup. You can find the software needed for download and execution here. For first-time users, a demo is available in the `demo` folder.

## Requirements
- Python version 3.7 or higher
- Libraries Required:
  - numpy
  - gzip
  - tqdm
    
(Install these libraries before the mapping process)

## Method 1: Mapping Reads to Generate a cell-by-cPeak Matrix (.mtx)

This method involves mapping sequencing reads from the `fragments.tsv.gz` file of each sample or cell to create a cell-by-cPeak matrix in .mtx format.

**Important**: Ensure that at least 20GB of memory is allocated to run this process using `--mode 'normal'`. For a quicker execution, switch to the 'performance' mode if you have at least 30GB of memory.

### Usage Instructions:

#### Step 1
Navigate to the `map2cpeaks` directory:
```bash
cd map2cpeaks
```
#### Step 2 To map your sequencing reads:

```bash
python main.py -f path/to/your_fragment.tsv.gz
```
• --fragment_path, -f: Input file path (must be a .gz file).

Optional Arguments:

• --help, -h: Display help information.
• --barcode_path, -b: Specify a file containing barcodes to be used. If not specified, all barcodes in the fragment file will be used.
• --output, -o: Designate an output folder (default is ./map2cpeaks_result).
• --output_name: Name the output files (default is cell-cpeaks).
• --mode: Choose between 'performance' (uses more memory for faster results) or 'normal' mode (default is 'performance').
• --reference: Specify the cPeak version, either 'hg38' or 'hg19' (default is 'hg38').

For example:

To run a demo mapping:

```bash
python main.py -f demo/test_fragment.tsv.gz
```
To use a provided barcode mapping (ensure 'barcodes.txt' is included in the fragments):

```bash
python main.py -f demo/test_fragment.tsv.gz -b demo/test_barcodes.txt --reference hg19
```
The resulting output will include a 'barcode.txt' and a '.mtx' file housing the mapping matrix.

To use hg19 as a reference of mapping, you can run:

```bash
python main.py -f demo/test_fragment.tsv.gz --reference hg19
```
## Method 2: Mapping Pre-Identified Features to cPeaks (Not Recommended)

Caution: This method can result in loss of genomic information as it only considers the pre-identified features. Moreover, the quantification of cPeaks may not be accurate for bulk ATAC-seq data.

Usage for Pre-Identified Feature Mapping:

```bash
python main.py [--bed_path feature.bed]
```
• --bed_path, -bed: Input the .bed file of features, such as 'MACS2calledPeaks.bed'.

Remember to priorly adjust your operational environment according to the system requirements and ensure you’ve properly understood the process to achieve optimal outcomes.
