# </a>Contents

* [Updates in 2020.02.28](#S-what_is_new)

* [Introduction](#S-introduction)

* [Installation](#S-installation)

* [Quick Start](#S-quick-start)

* [Input Format](#S-input-format)

* [Program Descriptions](#S-program-description)

* [Citation](#S-citation)

* [Contact](#S-contact)

* [Update Information](#S-update)

## <a name="S-what_is_new"></a> Updates of 2020.02.28

* Removed hard limit on maximum sequence length.  The hard limit always causes segmentation fault when performing ultra long sequence alignments.

* Merge `mecat2pw` and `mecat2ref` into one single mapping tool `mecat2map`. `mecat2map` uses much less memory compared to `mecat2pw`.

* The candidate partition process now supports multiple CPU threads.

* Support multiple input forms (see [Input Format](#S-input-format)). 

# <a name="S-introduction"></a>Introduction

MECAT2 is an improved version of  [MECAT](https://github.com/xiaochuanle/MECAT). It is an ultra-fast and accurate Mapping, Error Correction and de novo Assembly Tools for single molecula sequencing (SMRT) reads.


MECAT2 consists of three modules:

* `mecat2map`, a fast and accurate alignment tool for SMRT reads

* `mecat2cns`, correct noisy reads based on their pairwise overlaps

* `fsa`, a [string graph](https://academic.oup.com/bioinformatics/article/21/suppl_2/ii79/227189) based assembly tool.


MECAT2 is written in C, C++, and perl. It is open source and distributed under the [GPLv3](http://www.gnu.org/licenses/gpl-3.0.html) license.

**Please note that MECAT2 no longer supports Nanopore raw reads. We have developed a Mapping, Error Correction and de novo Assembly Pipeline specifically for Nanopore Raw Reads `NECAT`. [follow this link to NECAT](https://github.com/xiaochuanle/NECAT).**



# <a name="S-installation"></a>Installation

We have tested `MECAT2` on CentOS release 7.3 and on Ubuntu 18.04.

* Step 1: Figure out where to install `MECAT2`.
We will install `MECAT2` and two other auxiliary tools `HDF5` and `dextract`. We first identify the directory in which we want to install them. As an example, I will install them in the directory `/home/chenying/smrt_asm`. So I first create this directory using the `mkdir` command and go to that directory:
(**The dollar sign `$` that preceeds the input is the promt printed by the shell.**)
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
``` shell
$ git clone https://github.com/xiaochuanle/MECAT2.git
$ cd MECAT2
$ make 
$ cd ..
```
After installation, all the executables are found in `${MECAT_PATH}/MECAT/Linux-amd64/bin`. The folder name `Linux-amd64` will vary in operating systems.

* Step 3: Add relative pathes
``` shell
$ export PATH=${MECAT_PATH}/MECAT/Linux-amd64/bin:$PATH
```


# <a name="S-quick-start"></a>Quick Start

**Before running `MECAT2`, don't forget to add binary paths to `PATH` (Step 3 of [Installation](#S-installation)).**


Here we take assemblying the genome of Ecoli as an example, to go through each step in order.  Details of each step are given in the next section.

* Step 1: Download dataset.
We download the raw reads  [ecoli_filtered.fastq.gz](http://gembox.cbcb.umd.edu/mhap/raw/ecoli_filtered.fastq.gz) into directory `${MECAT_PATH}/ecoli`
``` shell
$ mkdir -p ${MECAT_PATH}/ecoli
$ cd ${MECAT_PATH}/ecoli
$ wget http://gembox.cbcb.umd.edu/mhap/raw/ecoli_filtered.fastq.gz
```
After that, we get raw read file `${MECAT_PATH}/ecoli/ecoli_filtered.fastq.gz`:
``` shell
$ ls
ecoli_filtered.fastq.gz
```

* Step 2: Prepare config file
We create a config file template using the following command:
``` shell
$ mecat.pl config ecoli_config_file.txt
```
This command creates a config file `ecoli_config_file.txt`, which looks like
``` shell
PROJECT=
RAWREADS=
GENOME_SIZE=
THREADS=4
MIN_READ_LENGTH=2000
CNS_OVLP_OPTIONS="-kmer_size 13"
CNS_PCAN_OPTIONS="-p 100000 -k 100"
CNS_OPTIONS=""
CNS_OUTPUT_COVERAGE=30
TRIM_OVLP_OPTIONS="-skip_overhang"
TRIM_PM4_OPTIONS="-p 100000 -k 100"
TRIM_LCR_OPTIONS=""
TRIM_SR_OPTIONS=""
ASM_OVLP_OPTIONS=""
FSA_OL_FILTER_OPTIONS="--max_overhang=-1 --min_identity=-1"
FSA_ASSEMBLE_OPTIONS=""
CLEANUP=0
```
After filling the relative information, we have
``` shell
PROJECT=ecoli
RAWREADS=/home/chenying/smrt_asm/ecoli/ecoli_filtered.fastq
GENOME_SIZE=4800000
THREADS=4
MIN_READ_LENGTH=2000
CNS_OVLP_OPTIONS="-kmer_size 13"
CNS_PCAN_OPTIONS="-p 100000 -k 100"
CNS_OPTIONS=""
CNS_OUTPUT_COVERAGE=30
TRIM_OVLP_OPTIONS="-skip_overhang"
TRIM_PM4_OPTIONS="-p 100000 -k 100"
TRIM_LCR_OPTIONS=""
TRIM_SR_OPTIONS=""
ASM_OVLP_OPTIONS=""
FSA_OL_FILTER_OPTIONS="--max_overhang=-1 --min_identity=-1"
FSA_ASSEMBLE_OPTIONS=""
CLEANUP=0
```

* Step 3: Correct Raw Reads.
Correct the raw noisy reads using the following command:
``` shell
$ mecat.pl correct ecoli_config_file.txt
```

* Step 4: Trim Out Low Quality Subsequences in Corrected Reads.
``` shell
$ mecat.pl trim ecoli_config_file.txt
```

* Step 5: Assemble Contigs Using the Trimeed Reads
``` shell
$ mecat.pl assemble ecoli_config_file.txt
```

* Step 6: Where to Find Results

	* The file `${MECAT_PATH}/eocli/ecoli/1-consensus/cns_reads_list.txt` contains the full path of all corrected reads files.
	``` shell
    $ cat ${MECAT_PATH}/eocli/ecoli/1-consensus//cns_reads_list.txt
    /home/chenying/smrt_asm/ecoli/ecoli/1-consensus/cns_cns_dir/p00000000.cns.fasta
    ```
	* The extracted longest 30x (The number 30 is indidated by the `CNS_OUTPUT_COVERAGE` option in the config file) corrected reads used for trimming is `${MECAT_PATH}/ecoli/ecoli/1-consensus/cns_final.fasta`.
	* The trimmed reads is `${MECAT_PATH}/ecoli/ecoli/2-trim_bases/trimReads.fasta`
	* The assembled contigs is `${MECAT_PATH}/ecoli/ecoli/4-fsa/contigs.fasta`



# <a name="S-input-format"></a>Input Format

The input to `MECAT2` is indicated by the `RAWREADS` option in the config file. It must be a **full path**. `MECAT2` supports several different input formats:

* `H5` format. `H5` file format must first be transfered to `FASTA` format with `${MECAT_PATH}/DEXTRACT/dextract`. For example:
```shell
$ find pathto/raw_reads -name "*.bax.h5" -exec readlink -f {} \; > reads.fofn
$ while read line; do   dextract -v $line >> reads.fasta ; done <  reads.fofn
```
After transformation, proceed to one of the following input case.

* `FASTA` format.
```shell
RAWREADS=/Users/sysu/Desktop/files/programs/ecoli/pacbio/ecoli/raw_reads.fasta
```
Or `FASTA` format compressed in [GNU Zip (gzip) format](https://www.gnu.org/software/gzip/manual/gzip.html)
```shell
RAWREADS=/Users/sysu/Desktop/files/programs/ecoli/pacbio/ecoli/raw_reads.fasta.gz
```

* `FASTQ` format
```shell
RAWREADS=/Users/sysu/Desktop/files/programs/ecoli/pacbio/ecoli/raw_reads.fastq
```
Or `FASTQ` format compressed in [GNU Zip (gzip) format](https://www.gnu.org/software/gzip/manual/gzip.html)
```shell
RAWREADS=/Users/sysu/Desktop/files/programs/ecoli/pacbio/ecoli/raw_reads.fastq.gz
```

* List format
A file indicates the **full paths** of all raw reads files.
``` shell
RAWREADS=/Users/sysu/Desktop/files/programs/tomato/read_list.txt
```
``` shell
$ cat /Users/sysu/Desktop/files/programs/tomato/read_list.txt
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161027_Spenn_001_001_all.fastq
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161101_Spenn_002_002_all.fastq
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161103_Spenn_003_003_all.fastq
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161108_Spenn_004_004_all.fastq
/share/home/chuanlex/xiaochuanle/data/testdata/tomato/20161108_Spenn_004_005_all.fastq
```
Please note that files in `read_list.txt` need not be the same format. Each file can independently be either `FASTA` or `FASTQ`, and can further be compressed in  [GNU Zip (gzip) format](https://www.gnu.org/software/gzip/manual/gzip.html).



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
CNS_OVLP_OPTIONS="-kmer_size 13"
CNS_PCAN_OPTIONS="-p 100000 -k 100"
CNS_OPTIONS=""
CNS_OUTPUT_COVERAGE=30
TRIM_OVLP_OPTIONS="-skip_overhang"
TRIM_PM4_OPTIONS="-p 100000 -k 100"
TRIM_LCR_OPTIONS=""
TRIM_SR_OPTIONS=""
ASM_OVLP_OPTIONS=""
FSA_OL_FILTER_OPTIONS="--max_overhang=-1 --min_identity=-1"
FSA_ASSEMBLE_OPTIONS=""
CLEANUP=0
```
The meaning of each option is given below

* `PROJECT=ecoli`, the name of the project. In this example, a directory `ecoli` will be created in the current directory, and then everything will take place in the directory `ecoli`.
* `RAWREADS=`, the raw reads (with full path) to be processed by `MECAT2`. See [Input Format](#S-input-format).
* `GENOME_SIZE=`, the size (in bp) of the underlying genome.
* `THREADS=`, number of CPU threads used  by `MECAT2`.
* `MIN_READ_LENGTH=`, minimal length of corrected reads and trimmed reads.
* `CNS_OVLP_OPTIONS=""`, options for detecting overlap candidates in the correction stage. Run `mecat2map -help` for details. Note that the output format is `seqidx` (`-outfmt seqidx`), which is set internally by `mecat.pl`. 
* `CNS_OPTIONS=""`, options for correcting raw reads. Run `mecat2cns -help` for details.
* `TRIM_OVLP_OPTIONS=""`, options for detecting overlaps in the trimming stage. Run `mecat2map` for details. Note that output format is `m4x` (`-outfmt m4x`), which is set internally by `mecat.pl`.
* `ASM_OVLP_OPTIONS=""`, options for detecting overlaps in the assemble stage. Run `mecat2map -help` for details. The output format is `m4` (`-outfmt m4`), which is set internally by `mecat.pl`.
* `FSA_OL_FILTER_OPTIONS=""`, options for filtering overlaps. See below for details.
* `FSA_ASSEMBLE_OPTIONS=""`, options for assembling trimmed reads. See below for details.
* `USE_GRID=false`, using multiple computing nodes (`true`) or not (`false`).
* `CLEANUP=0`, delete intermediate date genrated by `MECAT2` (`1`) or not (`0`). Please note the in assemblying large genomes, the intermediate data can be very large.
* `CNS_OUTPUT_COVERAGE=30`, number of coverage of the longest corrected reads are extracted to be trimed and then assembled. In this example, 30x (specifically, 30 * 4800000 = 144 MB) of the longest corrected reads will be extracted.

## <a name="SS-mecat2-workflow"></a> The `MECAT2` Workflow.

For easy use, we have integrated all the procedures into one perl script file `mecat.pl`, which works in the following steps:

> `meat.pl config`, as mentioned above, this command creates a config file.
> `mecat.pl correct`, correct raw reads, which consits of three steps:
>> detecting overlap candidates using `mecat2map`.
>>  partition overlap candidates into several parts using `mecat2pcan`. Each parts contains overlap candidates needed for correcting 100000 raw reads.
>>  correct raw reads based on overlap candidates using `mecat2cns`.

> `mecat.pl assemble`, assemble corrected reads in three steps:
>> extract 30x longest corrected reads with `mecat2extseqs`
>> trim out low quality subsequences in two stpes:
>>> detecting overlaps of extracted reads using `mecat2map`
>>> trim out low quality subsequence based on their overlaps using `mecat2lcr`, `mecat2splitreads` and `mecat2trimbases`.

>> assemble trimmed reads into contigs in three steps:
>>> detecting overlaps of trimmed reads using `mecat2map`
>>> filter out low quality overlaps using `fsa_ol_filter`
>>> assemble trimmed reads into contigs based on high quality overlaps using `fsa_assemble`


## <a name="SS-pairwise-mapping"></a> SMRT sequence alignment tool `mecat2map`



### </a>options

The command for running `mecat2pw` is

```shell
mecat2map [OPTIONS] reads reference > results.m4
```


## <a name="SS-consensus"></a> Correction Tool `mecat2cns`

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

Updates in MECAT2 (20193.14):

* Add some improvements in FSA

* Optimize Install Method

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

















