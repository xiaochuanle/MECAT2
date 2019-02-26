# </a>Contents

* [Introduction](#S-introduction)

* [Installation](#S-installation)

* [Quick Start](#S-quick-start)

* [Input Format](#S-input-format)

* [Program Descriptions](#S-program-description)

* [Citation](#S-citation)

* [Contact](#S-contact)

* [Update Information](#S-update)

# <a name="S-introduction"></a>Introduction

MECAT2 is an improved version of  [MECAT](https://github.com/xiaochuanle/MECAT). It is an ultra-fast and accurate Mapping, Error Correction and de novo Assembly Tools for single molecula sequencing (SMRT) reads. MECAT2 powers following improvements over MECAT:

* All modules in MECAT2 are developed by ourselves for easy maintenance.
* MECAT2 can be parallel to multiple nodes.
* Many bugs in MECAT have been fixed.
* To improve the mapping  efficiency, the default length of kmers is set to 15 (13 in MECAT).
* The default parameters in `mecat2cns` are relexed so that MECAT2 also works well for low coverage datasets.
* A new assembly module `fsa`, which is based on [string graph](https://academic.oup.com/bioinformatics/article/21/suppl_2/ii79/227189), is used for replacing the original assembly tool `mecat2canu`.


MECAT2 consists of four modules:

* `mecat2pw`, a fast and accurate pairwise mapping tool for SMRT reads

* `mecat2ref`, a fast and accurate reference mapping tool for SMRT reads

* `mecat2cns`, correct noisy reads based on their pairwise overlaps

* `fsa`, a [string graph](https://academic.oup.com/bioinformatics/article/21/suppl_2/ii79/227189) based assembly tool.


MECAT2 is written in C, C++, and perl. It is open source and distributed under the [GPLv3](http://www.gnu.org/licenses/gpl-3.0.html) license.

**Please note that MECAT2 no longer supports Nanopore raw reads. We have developed a Mapping, Error Correction and de novo Assembly Pipeline specifically for Nanopore Raw Reads `NECAT`, which will be released very soon.**



# <a name="S-installation"></a>Installation

We have tested `MECAT2` on CentOS release 7.3 and on Ubuntu 18.04.

* Step 1: Figure out where to install `MECAT2`.
We will install `MECAT2` and two other auxiliary tools `HDF5` and `dextract`. We first identify the directory in which we want to install them. As an example, I will install them in the directory `/home/chenying/smrt_asm`. So I first create this directory using the `mkdir` command and go to that directory:
(**The dollar sign `$` that preceeds the input is the promt that printed by the shell.**)
``` shell
$ mkdir -p /home/chenying/smrt_asm
$ cd /home/chenying/smrt_asm
$ pwd
/home/chenying/smrt_asm
```

For easy reference, we asign `/home/chenying/smrt_asm` to an environment variable `MECAT_PATH`:
``` shell
$ export MECAT_PATH=/home/chenying/smrt_asm
$ echo ${MECAT_PATH}
/home/chenying/smrt_asm
```

* Step 2:  Install `MECAT2`:

There are two ways to install `MECAT2`. One is to download executable binaries directly:
```shell
wget https://github.com/xiaochuanle/MECAT2/releases/download/20192026/mecat2_20190226_linuax_amd64.tar.gz
tar xzvf mecat2_20190226_linuax_amd64.tar.gz
```
Another way is to build from source codes:
``` shell
$ git clone https://github.com/xiaochuanle/MECAT2.git
$ cd MECAT2
$ make 
$ cd ..
```
After installation, all the executables are found in `${MECAT_PATH}/MECAT/Linux-amd64/bin`. The folder name `Linux-amd64` will vary in operating systems.

* Step 3: Install `HDF5`:
``` shell
$ wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.15-patch1/src/hdf5-1.8.15-patch1.tar.gz
$ tar xzvf hdf5-1.8.15-patch1.tar.gz
$ mkdir hdf5
$ cd hdf5-1.8.15-patch1
$ ./configure --enable-cxx --prefix=${MECAT_PATH}/hdf5
$ make
$ make install
$ cd ..
```
The header files of `HDF5` are put in `${MECAT_PATH}/hdf5/include`. The library files of `HDF5` are put in `${MECAT_PATH}/hdf5/lib` (in some systems, they are put in `${MECAT_PATH}/hdf5/lib64`, check it!).

* Step 4:  Install `dextract`
``` shell
$ cd DEXTRACTOR
$ export HDF5_INCLUDE=${MECAT_PATH}/hdf5/include
$ export HDF5_LIB=${MECAT_PATH}/hdf5/lib
$ make
$ cd ..
```

* Step 5: Add relative pathes
``` shell
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MECAT_PATH}/hdf5/lib
$ export PATH=${MECAT_PATH}/MECAT/Linux-amd64/bin:$PATH
$ export PATH=${MECAT_PATH}/DEXTRACTOR:$PATH
```


# <a name="S-quick-start"></a>Quick Start

**Before running `MECAT2`, don't forget to add library and binary paths to `LD_LIBRARY_PATH` and `PATH` (Step 5 of [Installation](#S-installation)).**


Here we take assemblying the genome of Ecoli as an example, to go through each step in order.  Details of each step are given in the next section.

* Step 1: Download dataset.

We download the raw reads  [ecoli_filtered.fastq.gz](http://gembox.cbcb.umd.edu/mhap/raw/ecoli_filtered.fastq.gz) into directory `${MECAT_PATH}/ecoli` and decompress it:

``` shell
$ mkdir -p ${MECAT_PATH}/ecoli
$ cd ${MECAT_PATH}/ecoli
$ wget http://gembox.cbcb.umd.edu/mhap/raw/ecoli_filtered.fastq.gz
$ gunzip ecoli_filtered.fastq.gz
```
After decompression, we get raw read file `${MECAT_PATH}/ecoli/ecoli_filtered.fastq`:
``` shell
$ ls
ecoli_filtered.fastq
```

* Step 2: Prepare config file

We create a config file template using the following command:
``` shell
$ mecat.pl config ecoli_config_file.txt
```
The above command creates a config file `ecoli_config_file.txt`, which looks like
``` shell
PROJECT=
RAWREADS=
GENOME_SIZE=
THREADS=4
MIN_READ_LENGTH=500
CNS_OVLP_OPTIONS=""
CNS_OPTIONS="-r 0.6 -a 1000 -c 4 -l 2000"
TRIM_OVLP_OPTIONS="-B"
ASM_OVLP_OPTIONS="-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400"
FSA_OL_FILTER_OPTIONS="--max_overhang=-1 --min_identity=-1"
FSA_ASSEMBLE_OPTIONS=""
USE_GRID=false
CLEANUP=0
CNS_OUTPUT_COVERAGE=30
GRID_NODE=0
```
After filling the relative information, we have

``` shell
PROJECT=ecoli
RAWREADS=/home/chenying/smrt_asm/ecoli/ecoli_filtered.fastq
GENOME_SIZE=4800000
THREADS=4
MIN_READ_LENGTH=2000
CNS_OVLP_OPTIONS=""
CNS_OPTIONS="-r 0.6 -a 1000 -c 4 -l 2000"
TRIM_OVLP_OPTIONS="-B"
ASM_OVLP_OPTIONS="-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400"
FSA_OL_FILTER_OPTIONS="--max_overhang=-1 --min_identity=-1"
FSA_ASSEMBLE_OPTIONS=""
USE_GRID=false
CLEANUP=0
CNS_OUTPUT_COVERAGE=30
GRID_NODE=0
```

* Step 3: Correct Raw Reads.

Correct the raw noisy reads using the following command:

``` shell
$ mecat.pl correct ecoli_config_file.txt
```

* Step 4: Assemble Contigs Using the Corrected Reads

``` shell
$ mecat.pl assemble ecoli_config_file.txt
```

* Step 5: Where to Find Results

	* The corrected reads is `${MECAT_PATH}/ecoli/ecoli/1-consensus/cns_reads.fasta`.
	* The extracted longest 30x corrected reads used for trimming is `${MECAT_PATH}/ecoli/ecoli/1-consensus/cns_final.fasta`.
	* The trimmed reads is `${MECAT_PATH}/ecoli/ecoli/2-trim_bases/trimReads.fasta`
	* The assembled contigs is `${MECAT_PATH}/ecoli/ecoli/4-fsa/contigs.fasta`



# <a name="S-input-format"></a>Input Format

`MECAT2` is capable of processing `FASTA`, `FASTQ`, and `H5` format files. However, the `H5` files must first be transfered to `FASTA` format by running `${MECAT_PATH}/DEXTRACT/dextract` before running `MECAT2`. For example:
```shell
$ find pathto/raw_reads -name "*.bax.h5" -exec readlink -f {} \; > reads.fofn
$ while read line; do   dextract -v $line >> reads.fasta ; done <  reads.fofn
```


# <a name="S-program-description"></a>Program Descriptions



We describe in detail each module of MECAT, including their options and output formats.


## <a name="SS-config-file"></a> Config File

`MECAT2` reads all the information, including project name, raw reads, and various running parameters, from config file. To create a config file template, just run

``` shell
$ mecat.pl config config_file_name
```
The above command creates a config file named `config_file_name`. We have met an sample of config file in the previous section

``` shell
PROJECT=ecoli
RAWREADS=/home/chenying/smrt_asm/ecoli/ecoli_filtered.fastq
GENOME_SIZE=4800000
THREADS=4
MIN_READ_LENGTH=2000
CNS_OVLP_OPTIONS=""
CNS_OPTIONS="-r 0.6 -a 1000 -c 4 -l 2000"
TRIM_OVLP_OPTIONS="-B"
ASM_OVLP_OPTIONS="-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400"
FSA_OL_FILTER_OPTIONS="--max_overhang=-1 --min_identity=-1"
FSA_ASSEMBLE_OPTIONS=""
USE_GRID=false
CLEANUP=0
CNS_OUTPUT_COVERAGE=30
GRID_NODE=0
```
The meaning of each option is given below

* `PROJECT=ecoli`, the name of the project. In this example, a directory `ecoli` will be created in the current directory, and then everything will take place in the directory `ecoli`.
* `RAWREADS=`, the raw reads (with full path) to be processed by `MECAT2`. Raw reads must be either FASTA or FASTQ format. The H5 format file must be transferred to FASTA format first (see [Input Format](#S-input-format)).
* `GENOME_SIZE=`, the size (in bp) of the underlying genome.
* `THREADS=`, number of CPU threads used  by `MECAT2`.
* `MIN_READ_LENGTH=`, minimal length of corrected reads and trimmed reads.
* `CNS_OVLP_OPTIONS=""`, options for detecting overlap candidates in the correction stage.
* `CNS_OPTIONS=""`, options for correcting raw reads.
* `TRIM_OVLP_OPTIONS=""`, options for detecting overlaps in the trimming stage.
* `ASM_OVLP_OPTIONS=""`, options for detecting overlaps in the assemble stage.
* `FSA_OL_FILTER_OPTIONS=""`, options for filtering overlaps.
* `FSA_ASSEMBLE_OPTIONS=""`, options for assembling trimmed reads.
* `USE_GRID=false`, using multiple computing nodes (`true`) or not (`false`).
* `CLEANUP=0`, delete intermediate date genrated by `MECAT2` (`1`) or not (`0`). Please note the in assemblying large genomes, the intermediate data can be very large.
* `CNS_OUTPUT_COVERAGE=30`, number of coverage of the longest corrected reads are extracted to be trimed and then assembled. In this example, 30x (specifically, 30 * 4800000 = 144 MB) of the longest corrected reads will be extracted.
* `GRID_NODE=0`, number of computing nodes used. This options is used only when `USE_GRID=1`.

## <a name="SS-mecat2-workflow"></a> The `MECAT2` Workflow.

For easy use, we have integrated all the procedures into one perl script file `mecat.pl`, which works in the following steps:

> `meat.pl config`, as mentioned above, this command creates a config file.
> `mecat.pl correct`, correct raw reads, which consits of two steps:
>> detecting overlap candidates using `mecat2pw`.
>>  correct raw reads based on overlap candidates using `mecat2cns`.

> `mecat.pl assemble`, assemble corrected reads in three steps:
>> extract 30x longest corrected reads
>> trim out low quality subsequences using `v2trim.sh` in two stpes:
>>> detecting overlaps of extracted reads using `v2asmpm.sh`
>>> trim out low quality subsequence based on their overlaps using `v2lcr` and `v2sr`.

>> assemble trimmed reads into contigs in three steps:
>>> detecting overlaps of trimmed reads using `v2asmpm.sh`
>>> filter out low quality overlaps using `fsa_ol_filter`
>>> assemble trimmed reads into contigs based on high quality overlaps using `fsa_assemble`


## <a name="SS-pairwise-mapping"></a> Pairwise Mapping Tool `mecat2pw`



### </a>options

The command for running `mecat2pw` is

```shell
mecat2pw -j [task] -d [fasta/fastq] -w [working folder] -t [# of threads] -o [output] -n [# of candidates] -a [overlap size] -k [# of kmers] -g [0/1]
```

The options are:

* `-j [task]`, job name, 0 = detect overlapping candidates only, 1 = output overlaps in M4 format, default = 1. If we are to correct noisy reads, outputing overlapping candidates is enough.

* `-d [fasta/fastq]`, reads file name in FASTA or FASTQ format.

* `-w [working folder]`, a directory for storing temporary results, will be created if not exists.

* `-t [# of threads]`, number of CPU threads used for overlapping, default=1.

* `-o [output]`, output file name

* `-n [# of candidates]`, number of candidates considered for gapped extension, default=100. Since each chunk  is about 2GB size, number of candidates(NC) should be set by genome size (GS).For GS < 20M, NC should be set as 200; For GS>20M and GS<200M; NC should be set as 100; For GS>200M, NC should be set as 50.

* `-a [overlap size]`, only output overlaps with length >= a. Default: 2000.

* `-k [# of kmers]`, two blocks between two reads having >= k kmer matches will be considered as a matched block pair. Default: 4.

* `-g [0/1]`, output the gapped extension start point (1) or not (0), default=0.


### </a>output format



If the job is detecting overlapping candidates, the results are output in `can` format, each result of which occupies one line and 9 fields:

```shell

[A ID] [B ID] [A strand] [B strand] [A gapped start] [B gapped start] [voting score] [A length] [B length]

```



`mecat2pw` outputs overlapping results in `M4` format, of which one result is given in one line. The fileds of `M4` format is given in order below:

```shell

[A ID] [B ID] [% identity] [voting score] [A strand] [A start] [A end] [A length] [B strand] [B start] [B end] [B length]

```

If the `-g` option is set to `1`, two more fields indicating the extension starting points are given:

```shell

[A ID] [B ID] [% identity] [voting score] [A strand] [A start] [A end] [A length] [B strand] [B start] [B end] [B length] [A ext start] [B ext start]

```

In the `strand` field, `0` stands for the forward strand and `1` stands for the reverse strand. All the positions are zero-based and are based on the forward strand, whatever which strand the sequence is mapped. Here are some examples:

```shell

44 500 83.6617 30 0 349 8616 24525 0 1 10081 21813

353 500 83.2585 28 0 10273 18410 22390 1 0 10025 21813

271 500 80.4192 13 0 14308 19585 22770 1 4547 10281 21813

327 501 89.8652 117 0 10002 22529 22529 1 9403 21810 21811

328 501 90.8777 93 0 0 10945 22521 1 0 10902 21811

```

In the examples above, read 500 overlaps with reads 44, 353, 271, 327 and 328.



### </a> memory consumption



Before overlapping is conducted, the reads will be split into several chunks. Each chunk is about 2GB in size so that the overlapping can be run on a 8GB RAM computer.



## <a name="SS-reference-mapping"></a> Reference Mapping Tool `mecat2ref`



### </a>options



`mecat2ref` is used for mapping SMRT reads to the reference genomes. The command  is

```shell

mecat2ref -d [reads] -r [reference] -w [folder] -t [# of threads] -o [output] -b [# of results] -m [output format]

```

The meanings of each option are as follows:

* `-d [reads]`, reads file name in FASTA/FASTQ format

* `-r [reference]`, reference genome file name in FASTA format

* `-w [folder]`, a directory for storing temporary results

* `-t [# of threads]`, number of working CPU threads

* `-o [output]`, output file name

* `-b [# of result]`, output the best b alignments

* `-m [output format]`, output format: 0 = ref, 1 = M4, 2 = SAM, default = 0


### </a>output format



`mecat2ref` outputs results in one of the three formats: the `ref` format, the `M4` format, and the `SAM` format.



For the `ref` format, each result occupies three lines in the form:

```shell

[read name] [ref name] [ref strand] [voting score] [read start] [read end] [read length] [ref start] [ref end]

mapped read subsequence

mapped reference subsequence

```

The strands of the reads are always forward. In the `[ref strand]` field, `F` indicates forward strand while `R` indicates reverse strand. All the positions are zero-based and relative to the forward strand. Here is an example:

```shell

1	gi|556503834|ref|NC_000913.3|	F	10	2	58	1988134	1988197

AAT-AGCGCCTGCCAGGCG-TCTTTT--CCGGCCATTGT-CGCAG--CACTGTAACGCGTAAAA

AATTAGCGCCTGCCAGGCGGTCTTTTTTCCGGCCATTGTTCGCAGGG-ACTGTAACGCGTAAAA

```

In this example, read 1 is mapped to the reference `gi|556503834|ref|NC_000913.3|`.



### </a>memory consumption



* Index for the genome: genomeSize * 8 bytes

* Compressed index for each CPU thread: genomeSize \* 0.1 * t bytes

* Local alignment: 100M * t + 1G bytes



## <a name="SS-consensus"></a> Correction Tool `mecat2cns`



`mecat2cns` is an adaptive error correction tool for high-noise single-molecula sequencing reads. It is as accurate as [pbdagcon](https://github.com/PacificBiosciences/pbdagcon.git) and as fast as FalconSense. Inputs to `mecat2cns` can be either `can` format or `M4` format. The command for running `mecat2cns` is

```shell

mecat2cns [options] overlaps-file reads output

```

The options are

* `-i [input type]`, input format, 0 = `can`, 1 = `M4`

* `-t [# of threads]`, number of CPU threads for consensus

* `-p [batch size]`, batch size the reads will be partitioned

* `-r [ratio]`, minimum mapping ratio

* `-a [overlap size]`, overlaps with length >= a will be used.

* `-c [coverage]`, minimum coverage, default=6

* `-l [length]`, minimum length of the corrected sequence

The default values for the other options are:
```shell
-i 1 -t 1 -p 100000 -r 0.6 -a 1000 -c 4 -l 2000
```


If the inputs are `M4` format, the overlap results in `[overlaps-file]` must contain the gapped extension start point, which means the option `-g` in `mecat2pw` must be set to 1, otherwise `mecat2cns` will fail to run. Also note that the memory requirement of `mecat2cns` is about 1/4 of the total size of the reads. For example, if the reads are of total size 1GB, then `mecat2cns` will occupy about 250MB memory.



### </a> output format



The corrected sequences are given in FASTA format. The header of each corrected sequence consists of three components seperated by underlines:

```shell

>A_B_C_D

```

where

* `A` is the original read id

* `B` is the left-most effective position

* `C` is the right-most effective position

* `D` is the length of the corrected sequence



by effective position we mean the position in the original sequence that is covered by at least `c` (the argument to the option `-c`) reads.

## <a name="SS-Extract longest correted reads"></a> `mecat2elr`



`mecat2elr` is used by `mecat.pl` for  extracting 30X longest sequences from the corrected data. The command is

``` shell

mecat2elr [the input fasta file from mecat2cns] [genome size]  [coverage] [the output filename]

```

## <a name="SS_fsa_ov_filter"></a> Overlap Filter `fsa_ol_filter`
`fsa_ol_filter` is used for filtering out low-quality overlaps. The usage of `fsa_ol_filter1` is

```shell
fsa_ol_filter [optioins] overlaps filtered_overlaps
```

The options are

* `--min_length=INT`, minimum length of reads (default: 2500)
* `--max_length=INT`, maximum length of reads (defualt: INT_MAX).
* `--min_identity=DOUBLE`, minimum identity of overlaps (defualt: 90).  
* `--min_aligned_length=INT`, minimum aligned length of overlaps (default: 2500). 
* `--max_overhang=INT`, maximum overhang of overlaps (default: 10), negative number = determined by the program.
* `--min_coverage=INT`, minimum base coverage (default: -1), negative number = determined by the program.
* `--max_coverage=INT`, maximum base coverage (default: -1), negative number = determined by the program.
* `--max_diff_coverage=INT`, maximum difference of base coverage (default: -1), negative number = determined by the program.
* `--coverage_discard=DOUBLE`, discard ratio of base coverage (default: 0.01). If `--max_coverage` or `--max_diff_coverage` is negative, it will be reset to (100-`coverage_discard`)th percentile.
* `--overlap_file_type="|m4|paf|ovl"`,  overlap file format (default: ""). `""` = filename extension, `"m4"` = `M4` format, `"paf"` = `PAF` format generated by minimap2, `"ovl"` = `OVL` format generated by FALCON.
* `--bestn=INT`, output best n overlaps on 5' or 3' end for each read (default: 10).
* `--genome_size=INT`, genome size. It determines the maximum length of reads with `--coverage` together. 
* `--coverage=INT`, coverage. It determines the maximum length of reads with `--genome_size` together.
* `--output_directory=STRING`, directory for output files (default: ".").
* `--thread_size=INT`, number of threads (default: 4).
           

## <a name="SS_fsa_assmble"></a> Assemble Tool `fsa_assemble`
`fsa_assemble` is a tool for constructing contigs from filtered overlaps and corrected reads. The algorithm is similar to FALCON. The usage of `fsa_assemble` is


```shell
fsa_assenble [optioins] filtered_overlaps
```

The options are
* `--min_length=INT`, minimum length of reads (default: 0).
* `--min_identity=DOUBLE`, minimum identity of overlaps (defualt: 0).
* `--min_aligned_length=INT`, minimum aligned length of overlaps (default: 0).
* `--min_contig_length=INT`, minimum length of contigs (default: 500). 
* `--read_file=STRING`, reads file name in FASTA or FASTQ format.
* `--overlap_file_type="|m4|paf|ovl"`,  overlap file format (default: ""). `""` = filename extension, `"m4"` = `M4` format, `"paf"` = `PAF` format generated by minimap2, `"ovl"` = `OVL` format generated by FALCON.
* `--output_directory=STRING`, directory for output files (default: ".").
* `--select_branch="no|best"`, selecting method when encountering branches in the graph, `"no"` = do not select any branch, `"best"` = select the most probable branch.     
* `--thread_size=INT`, number of threads (default: 4)
           

# <a name="S-citation"></a>`Citation`

Chuan-Le Xiao, Ying Chen, Shang-Qian Xie, Kai-Ning Chen, Yan Wang, Yue Han, Feng Luo, Zhi Xie. [MECAT: fast mapping, error correction, and de novo assembly for single-molecule sequencing reads](https://www.nature.com/articles/nmeth.4432). Nature Methods, 2017, 14: 1072-1074

# <a name="S-contact"></a> `Contact`

* Chuan-Le Xiao, xiaochuanle@126.com

* Ying Chen, chenying2016@gmail.com

* Fan Nie, niefan@csu.edu.cn

* Feng Luo, luofeng@clemson.edu


# <a name="S-update"></a>`Update Information`

Updates in MECAT2 (2019.2):

* Fix many bugs in MECAT

* Replace the asseble module mecat2canu by fasa.

Updates in MECAT V1.3 (2017.12.18):

* Correct text error in HDF5 Installation.

* Update the makefile in dextract .

* Update citation.

Updates in MECAT V1.2 (2017.5.22):

* Add `trimming module` in `mecat2canu` to improve the integrality of the assembly.

* Add supports for Nanopore data.

* Improve the sensitivity of `mecat2ref`.

MECAT v1.1 replaced the old MECAT,some debug were resolved and some new fuctions were added:
* 1. we added the extracted tools for the raw `H5` format files.
* 2. some debugs from running mecat2canu were solved

















