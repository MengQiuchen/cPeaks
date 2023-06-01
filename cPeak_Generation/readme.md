This function is uesd to generate the peak distribution from peak.bed

data: 2023/06/01


## Input 
The input of the code is bed files path, and you could use mclapply to calculate parallelly.

## output 

The output of the code is a list with 24 elements, each element is a dataframe. This dataframe contains "site" and "value" columns, "site" is the position of the peak region, and the "value" is always 1 because we only consider hitted positions in one bed.

## usage 

** Required: **

"path" : the bed to be transfered, it could contained many  bed files;
"savepath": the folder that you want to save the trasfered peak distribution

Optional: 
mc.cores: the number of threads.
chrlist: the focused chromsomes, default is chr 1:22, chrX and chrY.


## optional change
1. idx2list: a bin to a list, start-end,but **end -1**
2. chr2site: a chromsome range to dataframe with site and value 
3. file2site: a peak.bed to dataframe with site and value site

Each of the function could be changed if you want.

