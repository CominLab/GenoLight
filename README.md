# GenoLight
## Table of Contents

- [Getting Started](#started)
- [Users' Guide](#uguide)
  - [Installation](#install)
  - [General usage](#general)
  - [Getting help](#help)
[Citing GenoLight](#cite)


## <a name="uguide"></a>Users' Guide
Sequencing technologies have provided the basis of most modern genome sequencing studies due to their high base- level accuracy and relatively low cost. One of the central obstacles is mapping reads to the human reference genome. The reliance on a single reference human genome could introduce substantial biases in downstream analyses. More- over, including known variants in the reference makes read mapping, variant calling, and genotyping variant- aware. 
However, reads mapping is becoming a computationally intensive step for most genomic studies. Alignment-free methods have been used to save compute time and memory by avoiding the cost of full-scale alignment (Vinga and Almeida, Bioinformatics 2003). Recently, alignment-free approaches have been applied to SNP genotyping by (Shajiiet al., Bioinformatics 2016) and (Sun and Medvedev, Bioinformatics 2019). They introduce two SNP genotyping toolsnamed LAVA and 
VarGeno, respectively, which build an index from known SNPs (e.g. dbSNP) and then use approx-imate k-mer matching to genotype the donor from sequencing data. LAVA and VarGeno are reported to perform 4 to30 times faster than a standard alignment-based genotyping pipeline while achieving comparable accuracy. However,they require a large amount of memory, about 60GB. In this work, we introduce GenoLight that will address thisproblem with new efficient data structures.

---

## <a name="install"></a>Installation

To download and compile the code run the following commands.

First clone the repository, then cd into it and finally make it.
```sh
git clone --recursive https://github.com/frankandreace/GenoLight.git
cd GenoLight
make
```

You should now find the program in GenoLight folder.

---

##  <a name="general"></a>General Usage

GenoLight has two main functions: it creates a smart SNP dictionary from a reference genome and SNP dicitonary and then it performs the genotypization of a datasets.
The first phase is divided into 4 main functions:  

1. 'create_incomplete_smartSnpDictionary' needs 
    1. -r [reference genome]
    2. -s [snp dictionary]
    3. -p [output files prefix];

2. 'reassembly' needs
    1. -n [snp dictionary]
    2. -p [output files prefix];

3. 'createFMDIndex' doesn't need any additional command;

4. 'complete_smartSnpDictionary' needs
    1. -n [snp dictionary]
    2. -p [output files prefix].

Usage examples:

```sh
./genolight create_incomplete_smartSnpDictionary -r [reference genome] -s [snp dictionary] -p [output files prefix]
```
```sh
./genolight reassembly -n [snp dictionary] -p [output files prefix]
```
```sh
./genolight createFMDIndex   
```
```sh
./genolight complete_smartSnpDictionary -n [snp dictionary] -p [output files prefix]
```

The second one has only one function.

1.'geno' needs
    1.-t [number of threads]
    2.-r [single end reads]
    3.-rl [Left/.1 paired end reads]
    4.-rr [Right/.2 paired end reads]
    5.-q [quality score]
    6.-s [snp dictionary]
    7.-p [output files prefix]
    8.-o [genotyping output file].

NOTE THAT you should EITHER use '-r' or '-rl' & '-rr'.

Usage example: 

Paired end 

```sh
./genolight geno -t [number of threads] -rl [Left/.1 paired end reads] -rr [Right/.2 paired end reads] -s [snp dictionary] -p [output files prefix] -o [genotyping output file]
```
or 

Single end 

```sh
./genolight geno -t [number of threads] -r [single end reads] -s [snp dictionary] -p [output files prefix] -o [genotyping output file]
```
---

## <a name="help"></a>Getting help
If you encounter bugs or have further questions or
requests, you can raise an issue at the [issue page][issue]. You can also contact Francesco Andreace at francesco.andreace@unipd.it .

---

## <a name="cite"></a>Citing GenoLight

GenoLight paper has been accepted at BITS 2021.
If you use GenoLight in your work, please cite:

[issue]: https://github.com/frankandreace/GenoLight/issues
