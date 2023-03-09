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

**- Input files**

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






