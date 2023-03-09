# Allele_auto
## **A pipeline for allele identification and allele-specific gene expression.**

**Version:** 1.0.0

**Author:** Tian-Le Shi

**Email:** tianleshi@bjfu.edu.cn

**Link:** https://github.com/shitianle77/Allele_auto
***

## 1 What is Allele_auto

**Program:** Allele_auto: a pipeline for allele identification and allele-specific gene expression

**Main commands:**

**[ pipeline ]**

  allele_identification for allele identification
  
  allele_specific_expression for allele-specific expression analysis

**For details on usage, please initiate the following commands:**

  ```
  bash allele_identification.sh -h
  
  bash allele_specific_expression.sh -h
```

**[ Rule for allele-pair identification ]**

  3σ rule was used for allele-pair identification, please refer to its original publication for details [(Lehmann, 2013)](https://ascelibrary.org/doi/abs/10.1061/(ASCE)SU.1943-5428.0000112). 
  
## 2 Installation
### 2.1 Dependencies

Genetribe: https://github.com/chenym1/genetribe

WGDI: https://github.com/SunPengChuan/wgdi

#Please refer to the respective links, for the installation of these abovementioned softwares.

### 2.2 Installation

1. Download the latest version of Allele_auto from Github. We only provide Linux 64-bit version.

`git clone https://github.com/shitianle77/Allele_auto.git`

2. Adding the environment variables to your system PATH, if necessary.

Change your PATH following these steps.

`vi ~/.bash_profile`

add Genetribe, WGDI and R environments to ~/.bash_profile (eg: export PATH=/home/tlshi/software/genetribe:$PATH).

Then apply the changes in your PATH, with this command.

`source ~/.bash_profile`

Download and install the required R packages (DESeq2, ggplot2) by executing the following code in your R console.

`install.packages(c(“DESeq2”, “ggplot2”))`

3. Get helps on Allele_auto.
```
cd Allele_auto
bash ./bin/allele_identification.sh -h
bash ./bin/allele_specific_expression.sh -h
```

## 3 Quick Start

We provide example files in the test folder. Here we show steps on running a simple computation with this pipeline. We expect this may help users to quickly get familiar with this pipeline.

```
cd Allele_auto
ls
# 00_data 01_genetribe 02_wgdi 03_DEG bin
```

The software is divided into two main sections, allele identification and allele-specific gene expression. All the input data are under the directory 00_data. You can keep only the directories 00_data and bin and then execute the following lines to get the results of the example data.

Before the official start, users need to modify them:

(1)	Once genetribe has been successfully installed, you will need to go into the ~/genetribe/bin directory and modify the penultimate line in core to remove the “#”.

(2)	In the script ./bin/allele_identification.sh, 195 lines need to be replaced with the user's environment.

```
bash ./bin/allele_identification.sh -p chrpairs.txt -a SA -b SB
bash ./bin/allele_specific_expression.sh -a SA -b SB -c allelepairs.count_selected.txt -t allelepairs.tpm_selected.txt -s 102
```

## 4 Input and output

Here we have described in details on formats of the input files and the output files.

### 4.1 Files for the procedure of allele identification

- **Input files**

(1) One such input file is of the homologous chromosome grouping.

**Get an example of such input by executing this command:**

```
cat chrpairs.txt
# chr01A	chr01B	chr01
# chr02A	chr02B	chr02
```

(2) Other input files include those of protein sequences (for example: SA.pep, SB.pep), cds sequences (SA.cds, SB.cds), fasta sequences (SA.fa, SB.fa) and gff files (SA.gff3, SB.gff3) from gene annotation and genome assembly of both the two haplotype (two subgenomes for a haplotype-resolved genome assembly of a diploid species) genomes.

**Get examples of such input by executing the following commands:**

```
cat SA.pep
# >PAxG01Ag0000100
# MLGIIGLLETKVISAQMETVITNLQLPSWNFISNINASSSCRVIVGWDPTMFRISCLHNSDQ

cat SB.pep
# >PAxG01Bg0000100
# MHRAAGGFGGLRINSRLLIPVPYADPEDDYTVILNDWYTSSHATLRKYLDDGRSLARPTG

cat SA.cds
# >PAxG01Ag0000100
# ATGCTGGGTATCATTGGTCTGCTGGAAACTAAGGTTATATCAGCCCAAATGGAGACTGTTATCACTAATTTACAGCTACCTTCTTGGAATTTTATTTCAAATATAAATGCCTCCTCTAGCTGTCGGGTGATTGTGGGATGGGATCCTACAATGTTTCGTATTTCTTGTCTGCATAATTCT

cat SB.cds
# >PAxG01Bg0000100
# ATGCACCGTGCAGCTGGAGGCTTTGGCGGCCTTCGCATCAACAGCCGCCTACTCATCCCTGTACCTTATGCTGATCCCGAGGATGACTACACCGTCATACTTAATGACTGGTATACCAGCAGCCACGCCACTCTCAGGAAATACTTGGATGACGGCCGCTCTCTTGCAAGGCCTAC

cat SA.gff3
# chr01A  maker  gene  16571  20155  .  -  .  ID=PAxG01Ag0000100;Name=maker-chr01A:1-4000000-exonerate_protein2genome-gene-0.103;Alias=maker-chr01A:1-4000000
# chr01A  maker  mRNA  16571  20155  1782.0  -  .  ID=PAxG01Ag0000100.1;Parent=PAxG01Ag0000100;Name=maker-chr01A:1-4000000-exonerate_protein2genome-gene-0.103-mRNA-1
# chr01A  maker  exon  16571  17453  .  -  .  ID=maker-chr01A:1-4000000-exonerate_protein2genome-gene-0.103-mRNA-1:exon:2101;Parent=PAxG01Ag0000100.1
# chr01A  maker  exon  19257  20155  .  -  .  ID=maker-chr01A:1-4000000-exonerate_protein2genome-gene-0.103-mRNA-1:exon:2100;Parent=PAxG01Ag0000100.1
# chr01A  maker  CDS  19257  20155  .  -  0  ID=maker-chr01A:1-4000000-exonerate_protein2genome-gene-0.103-mRNA-1:cds;Parent=PAxG01Ag0000100.1
# chr01A  maker  CDS  16571  17453  .  -  1  ID=maker-chr01A:1-4000000-exonerate_protein2genome-gene-0.103-mRNA-1:cds;Parent=PAxG01Ag0000100.1

cat SB.gff3
# chr01B  maker  gene  39113  40364  .  -  .  ID=PAxG01Bg0000100;Name=maker-chr01B:1-4000000-augustus-gene-0.16;Alias=maker-chr01B:1-4000000
# chr01B  maker  mRNA  39113  40364  .  -  .  ID=PAxG01Bg0000100.1;Parent=PAxG01Bg0000100;Name=maker-chr01B:1-4000000-augustus-gene-0.16-mRNA-1
# chr01B  maker  exon  39522  40364  .  -  .  ID=maker-chr01B:1-4000000-augustus-gene-0.16-mRNA-1:exon:416;Parent=PAxG01Bg0000100.1
# chr01B  maker  exon  39113  39118  .  -  .  ID=maker-chr01B:1-4000000-augustus-gene-0.16-mRNA-1:exon:415;Parent=PAxG01Bg0000100.1
# chr01B  maker  CDS  39522  40364  .  -  0  ID=maker-chr01B:1-4000000-augustus-gene-0.16-mRNA-1:cds;Parent=PAxG01Bg0000100.1
# chr01B  maker  CDS  39113  39118  .  -  0  ID=maker-chr01B:1-4000000-augustus-gene-0.16-mRNA-1:cds;Parent=PAxG01Bg0000100.1

cat SA.fa
# >chr01A
# aaaccctaaaccctaaaccctaaaccctaaaCCCTAAACCCTAAACCCTAAACCCTAAaCCCTAAA…
# >chr02A
# acccTAAACCCTAAACCCTAAACCCTAAACCCTTAAACCCTAAACCCTAAACCCTA…

cat SB.fa
# >chr01B
# taaaccctaaaccctaaaccctaaaccctaaaccctaaaCCCTAAACCCTAAACCCTAAACCCTAAACC…
# >chr02B
# aaaccctaaaccctaaaccctaaaccctaaacccTAAACCCTAAACCCTAAACCCTAAACCCTAAAC… 
```
We also have to prepare input files of genes (in bed format, for example: SA.bed and SB.bed) and chromosome length information (SA.len and SB.len), which can be obtained from the genome assembly and annotation of the haplotype (two subgenomes for a haplotype-resolved genome assembly of a diploid species) genomes.

(1) prepare bed file from gff file of subgenomes A and B (This is done by using the JCVI)

```
python -m jcvi.formats.gff bed --type=gene --key=ID SA.gff3 -o SA.bed
python -m jcvi.formats.gff bed --type=gene --key=ID SB.gff3 -o SB.bed
```
Output files: SA.bed SB.bed. column are chromosome, start location, end location, gene ID, score, and strand.

**Note:** For some annotations, the fourth column of bed contains some strings, such as gene:, ID:, etc., which we should delete.

**Get an example of such input by executing the following commands:**

```
cat SA.bed
# chr01A	16570	20155	PAxG01Ag0000100	0	-
# chr01A	40601	47447	PAxG01Ag0000200	0	-

cat SB.bed
# chr01B	39112	40364	PAxG01Bg0000100	0	-
# chr01B	47342	49123	PAxG01Bg0000200	0	+
```

(2) prepare the file of chromosome length of two subgenomes A and B

```
python generate_conf.py -p SA SA.fa SA.gff3
python generate_conf.py -p SB SB.fa SB.gff3
```
Output files: SA.len SB.len SA.gff SB.gff (The gff files are used for subsequent filtering.)

```
cat SA.len
# chr01A	50794782	4038
# chr02A	24505410	2307

cat SB.len
# chr01B	51099886	4086
# chr02B	27466271	2296
```
Columns are chromosome, length of chromosome and number of chromosome genes.


```
cat SA.gff
# chr01A	PAxG01Ag0000100 16571	20155	-	1	PAxG01Ag0000100.1
# chr01A	PAxG01Ag0000200 40602	47447	-	2	PAxG01Ag0000200.1

cat SB.gff
# chr01B	PAxG01Bg0000100	39113	40364	-	1	PAxG01Bg0000100.1
# chr01B	PAxG01Bg0000200	47343	49123	+	2	PAxG01Bg0000200.1
```
Each column is chromosome number, gene name, start location, end location, strand, order of each chromosome (starting from 1) and original id and not read.

For more information, please refer to: https://wgdi.readthedocs.io/en/latest/usage.html 

- **Output files**

**The expected output files are written in different folders. Here, we have** “raw_RBH.genepairs”**,** “RBH.genepairs”**,** **and** “SA_SB.blast” **in the** “01_genetribe” **folder,** **and we have** “SA_SB.collinearity.txt”**,** “SA_SB.ks.txt”**,** “SA_SB_block.csv”**,** “block.tsv”**,** “filtered.block.tsv”**,** “genepairs_info.tsv”**,** “filtered.genepairs_info.tsv”**,** **and** “allele_pairs_lasted.txt” **in the** “02_wgdi” **folder.**

**Get into different folders for details, by executing the following commands:**
```
cd 01_genetribe
```
|Output files|Description|
|---|---|
|raw_RBH.genepairs|Gene pairs belonging to the Reciprocal Best Hits (RBH)|
|RBH.genepairs|List of gene pairs belonging to RBH|
|SA_SB.blast|Blast information for gene pairs belonging to RBH|

For more information, please refer to: https://chenym1.github.io/genetribe/tutorial/fileformats.html

```
cd 02_wgdi
```
|Output files|Description|
|---|---|
|SA_SB.collinearity.txt|Improved collinearity (For details, see https://wgdi.readthedocs.io/en/latest/collinearity.html)|
|SA_SB.ks.txt|Non-synonymous (Ka) and synonymous (Ks) (For details, see https://wgdi.readthedocs.io/en/latest/ks.html)|
|SA_SB_block.csv|BlockInfo (For details, see https://wgdi.readthedocs.io/en/latest/blockinfo.html)|
|block.tsv|BlockInfo (Same as SA_SB_block.csv, separated by tab)|
|filtered.block.tsv|Filtered blocks information|
|genepairs_info.tsv|Details of allele pairs on blocks (Position and Ks information)|
|filtered.genepairs_info.tsv|Details of filtered allele pairs on blocks (Position and Ks information)|
|allele_pairs_lasted.txt|The final allele pairs identified|

### 4.2 Files for computation of allele-specific expression

- **Input files**

The two core input files are of the Count and TPM expression matrices for the allele pairs.

The count matrix has the following format, and the TPM matrix has the same format.

|Output files|Description|
|---|---|
|SA_SB.collinearity.txt|Improved collinearity (For details, see https://wgdi.readthedocs.io/en/latest/collinearity.html)|


SA_s1_1: The first repeat (1) of the allele of subgenome A (SA) in sample1 (s1).

SB_s2_3: The third repeat (3) of the allele of subgenome B (SB) in sample2 (s2).

**Here is one example (The count values of allele pairs in the three replicates of the first sample.):**

```
cd ./00_data/RNA_seq
cat allelepairs.count_selected.txt
# PAxG01Ag0000300_PAxG01Bg0000100		0	0	0	0	0	0
# PAxG01Ag0000400_PAxG01Bg0000200		0	0	0	0	0	0
# PAxG01Ag0000700_PAxG01Bg0000400		780.324	906.261	796.327	712.676	832.739	733.67
# PAxG01Ag0000800_PAxG01Bg0000500		1533.223	1573.857	1557.604	898.777	1141.143	875.396
# PAxG01Ag0000900_PAxG01Bg0000600		0	0	0	0	14	1
```

- **Output files**
**The expected output files are written in different folders. Here, we have** “pairs.genelist”**,** “name_list.txt”**,** “AvsB.”**,** “stats_number/”**,** **and** “High_expression/” **in the** “03_DEG/1_Class_alleles” **folder, and we have** “.name”**,** “.tpm”**,** “tpm.box.txt”**,** “tpm_boxplot.pdf” **in the** “03_DEG/2_Diff_comparison” **folder, and also** “all.kaks” **and** “.pdf”**, in the “KaKs” folder.**


















