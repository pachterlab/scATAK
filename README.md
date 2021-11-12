## scATAK tool for pre-processing 10X single-cel ATAC-seq data.

### This repository includes
* 1) A BASH script scATAK that wraps several bioinformatics tools (e.g., kallisto, bustools) for pre-processing single-cell ATAC-seq data; 
* 2) A bin folder with Linux executables and a lib folder with the whitelists for 10X single-cell ATAC-seq data.

### Installation
* 1) Download scATAK
```
cd ~
git clone https://github.com/pachterlab/scATAK.git
```

* 2) Add `export SCATAK_HOME=/home/fgao/scATAK/` to `.bashrc` file, and then `source .bashrc`.  
