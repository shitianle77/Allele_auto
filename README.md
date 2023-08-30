# Allele_auto
## **A pipeline for allele identification and allelespecific gene expression with haplotype-resolved diploid genome assembly.**

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

(1) 3σ rule is used for allele-pair identification, please refer to its original publication for details [(Lehmann, 2013)](https://ascelibrary.org/doi/abs/10.1061/(ASCE)SU.1943-5428.0000112). 

(2) Tukey’s method is an acclaimed and straightforward graphical technique known forrepresenting continuous univariate data through a boxplot. This method calculates the upper and lower extremes of the data through the quartile.
  
## 2 Installation
### 2.1 Dependencies

Genetribe: https://github.com/chenym1/genetribe

WGDI: https://github.com/SunPengChuan/wgdi

#Please refer to the respective links, for the installation of these abovementioned software.

### 2.2 Installation

**2.2.1. Download the latest version of Allele_auto from Github. We only provide Linux 64-bit version.**

`# git clone https://github.com/shitianle77/Allele_auto.git`

**2.2.2. Adding the environment variables to your system PATH, if necessary.**
Change your PATH following these steps.

`# vi ~/.bash_profile`

add Genetribe, WGDI and R environments to ~/.bash_profile (eg: export PATH=/home/tls/genetribe:$PATH).
Then apply the changes in your PATH, with this command.

`# source ~/.bash_profile`

Download and install the required R packages (DESeq2, ggplot2) by executing the following code in your R console.

`# install.packages(c(“DESeq2”, “ggplot2”))`

**2.2.3. Get helps on Allele_auto.**
```
# cd Allele_auto
# bash ./bin/allele_identification.sh -h
# bash ./bin/allele_specific_expression.sh -h
```

## 3 Quick Start

We provide example files in the test folder. Here we show steps on running a simple computation with this
pipeline. We expect this may help users to quickly get familiar with this pipeline.

```
# cd Allele_auto
# ls
00_data 01_genetribe 02_wgdi 02_wgdi‐Tukey  03_DEG bin
```

The software is divided into two main sections, allele identification and allele-specificgene expression. All
the input data are under the directory 00_data. You can keep only the directories 00_data and bin and then
execute the following lines to get the results of the example data.

Before the official start, users need to modify them:

(1)	Once genetribe has been successfully installed, you will need to go into the ~/genetribe/bin directory and modify the penultimate line in core to add the “#”.
```
# sed -i 's/rm -rf output/#rm -rf output/'  yourgenetribepath/genetribe/bin/core
```
(2)	In the script ./bin/allele_identification.sh, 195 lines need to be replaced with the user's environment.

```
# bash ./bin/allele_identification.sh ‐p Zo_chrpairs.txt ‐a Zo_SA ‐b Zo_SB
# bash ./bin/allele_specific_expression.sh ‐a Zo_SA ‐b Zo_SB ‐c allelepairs.count_selected.txt ‐t allelepairs.tpm_selected.txt ‐s 21
```
**Note:** After obtaining the RBH gene pair (after running bash allele_identification.sh -p Zo_chrpairs.txt -a
Zo_SA -b Zo_SB), we also provided another way to filter alleles through a single line of code: bash
./bin/allele_filtering.sh - p Zo_chrpairs.txt - a Zo_SA- b Zo_SB.

## 4 Input and output

Here we have described in details on formats of the input files and the output files.

### 4.1 Files for the procedure of allele identification

- **Input files**

(1) One such input file is of the homologous chromosome grouping

**Get an example of such input by executing this command:**

```
# cat Zo_chrpairs.txt
chr01A	chr01B	chr01
chr02A	chr02B	chr02
```

(2) Other input files include those of protein sequences (for example: Zo_SA.pep, Zo_SB.pep), cds sequences (Zo_SA.cds, Zo_SB.cds), fasta sequences (Zo_SA.fa, Zo_SB.fa) and gff files (Zo_SA.gff3, Zo_SB.gff3) from gene annotation and genome assembly of both the two haplotype (two subgenomes for a haplotype-resolved genome assembly of a diploid species) genomes.

**Get examples of such input by executing the following commands:**

```
# cat Zo_SA.pep
>Zo_SA_chr01A_1
MNNASPSAAEPNSHALALPNPSSPLKDRSTYTNLKEHLLRPAGNNLWSPPVSKRATAGSK

# cat Zo_SB.pep
>Zo_SB_chr01B_1
MITTRFFPHSRFFLPSHLPTLCRPIHSGAAHPRITRSELVDRICRILTLERFHAIPKLPFRFSDD

# cat Zo_SA.cds
>Zo_SA_chr01A_1
ATGAATAATGCAAGTCCATCTGCTGCAGAACCCAACTCACACGCACTTGCTCTTCCTAATCCTTCTTCCCCACTTAAAGACCGTTCTACATATACTAACTTAAAAGAACATCTCTTGAGACCAGCAGGGAACAATCTCTGGTCGCCACCAGTCAGCAAGCGAGCGACAGCCGGGAG

# cat Zo_SB.cds
>Zo_SB_chr01B_1
TGATCACAACTAGATTCTTCCCTCATTCGCGTTTCTTCCTCCCCTCGCACCTGCCCACTCTCTGCCGGCCCATCCACTCCGGCGCTGCCCACCCCCGCATCACCAGATCTGAGCTCGTCGACCGGATATGCCGCATCCTCACCCTCGAGCGCTTCCACGCCATTCCCAAGCTTCCC

cat Zo_SA.gff3
# chr01A	maker	gene	1954	4056	.	+	.	ID=Zioff01G0000100;Name=Zioff01G0000100
# chr01A	maker	mRNA	1954	4056	.	+	.	ID=Zioff01G0000100.1;Parent=Zioff01G0000100;Name=Zioff01G0000100.1;
# chr01A	maker	exon	1954	2105	.	+	.	ID=Zioff01G0000100.1:exon:1170;Parent=Zioff01G0000100.1
# chr01A	maker	exon	3204	4056	.	+	.	ID=Zioff01G0000100.1:exon:1171;Parent=Zioff01G0000100.1
# chr01A	maker	CDS	1954	2105	.	+	0	ID=Zioff01G0000100.1:cds;Parent=Zioff01G0000100.1
# chr01A	maker	CDS	3204	4056	.	+	1	ID=Zioff01G0000100.1:cds;Parent=Zioff01G0000100.1

cat Zo_SB.gff3
# chr01B	maker	gene	2566	8399	.	+	.	ID=Zioff01G0466400;Name=Zioff01G0466400
# chr01B	maker	mRNA	2566	8399	.	+	.	ID=Zioff01G0466400.1;Parent=Zioff01G0466400;Name=Zioff01G0466400.1
# chr01B	maker	exon	2566	2954	.	+	.	ID=Zioff01G0466400.1:exon:1;Parent=Zioff01G0466400.1
# chr01B	maker	CDS	2680	2954	.	+	0	ID=Zioff01G0466400.1:cds;Parent=Zioff01G0466400.1
# chr01B	maker	CDS	3045	5487	.	+	1	ID=Zioff01G0466400.1:cds;Parent=Zioff01G0466400.1
# chr01B	maker	exon	3045	5487	.	+	.	ID=Zioff01G0466400.1:exon:2;Parent=Zioff01G0466400.1
# chr01B	maker	CDS	5845	5918	.	+	0	ID=Zioff01G0466400.1:cds;Parent=Zioff01G0466400.1
# chr01B	maker	exon	5845	5918	.	+	.	ID=Zioff01G0466400.1:exon:3;Parent=Zioff01G0466400.1
# chr01B	maker	CDS	7760	7949	.	+	1	ID=Zioff01G0466400.1:cds;Parent=Zioff01G0466400.1
# chr01B	maker	exon	7760	8399	.	+	.	ID=Zioff01G0466400.1:exon:4;Parent=Zioff01G0466400.1

cat Zo_SA.fa
# >chr01A
# taGCAAGTTGTTTTACCTAATTTATTTTAATGTTAAATATTTAGTATTTGTTGATAAAAA…
# >chr02A
# caccggcggccggaggagttggggcggattgcaggactttggacccagCGCCGGCGCCTCCTTCGCATGGGA…

cat Zo_SB.fa
# >chr01B
# ccaaatagttgatactacttgcccatgggtttcaaaggtatttgtttcccttttctaTCAGAGTAGAGAATAAGGTCTTG…
# >chr02B
# gtatgcttagtgtccagatgccaatgccggcggcaaggagtcatgggagcatctggctgggagcctcgcggtggagatgtccggggt…
```
We also have to prepare input files of genes (in bed format, for example: Zo_SA.bed and Zo_SB.bed) and chromosome length information (Zo_SA.len and Zo_SB.len), which can be obtained from the genome assembly and annotation of the haplotype (two subgenomes for a haplotype-resolved genome assembly of a diploid species) genomes.

(1) prepare bed file from gff file of subgenomes A and B (This is done by using the JCVI)

```
# python -m jcvi.formats.gff bed --type=gene --key=ID Zo_SA.gff3 -o Zo_SA.bed
# python -m jcvi.formats.gff bed --type=gene --key=ID Zo_SB.gff3 -o Zo_SB.bed
```
Output files: Zo_SA.bed and Zo_SB.bed. Columns are chromosome, start location, end location, gene ID, score, and strand.

**Note:** For some annotations, the fourth column of bed contains some strings, such as gene:, ID:, etc., which we should delete.

**Get an example of such input by executing the following commands:**

```
cat Zo_SA.bed
# chr01A	1953	4056	Zo_SA_chr01A_1	0	+
# chr01A	15432	24328	Zo_SA_chr01A_2	0	+

cat Zo_SB.bed
# chr01B	2565	8399	Zo_SB_chr01B_1	0	+
# chr01B	9304	13618	Zo_SB_chr01B_2	0	-
```

(2) prepare the file of chromosome length of two subgenomes A and B

```
python generate_conf.py -p Zo_SA Zo_SA.fa Zo_SA.gff3
python generate_conf.py -p Zo_SB Zo_SB.fa Zo_SB.gff3
```
Output files: Zo_SA.len Zo_SB.len Zo_SA.gff Zo_SB.gff (The gff files are used for subsequent filtering.)

```
cat Zo_SA.len
# chr01A	50794782	4038
# chr02A	24505410	2307

cat Zo_SB.len
# chr01B	171592823	4318
# chr02B	156697907	3810
```
Columns are chromosome, length of chromosome and number of chromosome genes.


```
cat Zo_SA.gff
# chr01A	Zo_SA_chr01A_1	1954	4056	+	1	Zo_SA_chr01A_1.1
# chr01A	Zo_SA_chr01A_2	15433	24328	+	2	Zo_SA_chr01A_2.1

cat Zo_SB.gff
# chr01B	Zo_SB_chr01B_1	2566	8399	+	1	Zo_SB_chr01B_1.1
# chr01B	Zo_SB_chr01B_2	9305	13618	-	2	Zo_SB_chr01B_2.1
```
##Each column is chromosome number, gene name, start location, end location, strand, order of each chromosome (starting from 1) and original id and not read.

##For more information, please refer to: https://wgdi.readthedocs.io/en/latest/usage.html 

**Note:** To facilitate the subsequent process, we have renamed the gene IDs, with the correspondence between the original and updated gene IDs provided in the two files (rename_SA.genelist and rename_SB.genelist). 

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

```
cd 02_wgdi-Tukey
```
|Output files|Description|
|---|---|
|allele_pairs_lasted.txt|The final allele pairs identified|
|block.tsv|BlockInfo (Same as SA_SB_block.csv, separated by tab)|
|genepairs_info.tsv|Details of allele pairs on blocks (Position and Ks information)|
|filtered_ks_NG86-1.5.genepairs_info.tsv|Details of filtered allele pairs after removing Ks outliers (Position and Ks information)|
|filtered_slope-1.5.genepairs_info.tsv|Details of filtered allele pairs after removing slope outliers (Position and Ks information)|
|filtered.ks-slopebioplot_genepairs_info.tsv|Details of filtered allele pairs after removing Ks and slope outliers (Position and Ks information)|

```
cd 02_wgdi-Tukey/allele_plot
```
|Output files|Description|
|---|---|
|chr01.coord.allele|Details of the position information of the filtered allele pairs on chromosome 1|
|chr02.coord.allele|Details of the position information of the filtered allele pairs on chromosome 2|
|pairs.coord.allele|Details of the position information of the filtered allele pairs|
|pairs.allele.pdf|Dotplot plot of filtered allele pairs|

### 4.2 Files for computation of allele-specific expression

- **Input files**

The two core input files are of the Count and TPM expression matrices for the allele pairs.

The count matrix has the following format, and the TPM matrix has the same format.

|Allele_ID|SA_s1_1|SA_s1_2|SA_s1_3|SA_s2_1|SA_s2_2|SA_s2_3|...|SB_s1_1|SB_s1_2|SB_s1_3|SB_s2_1|SB_s2_2|SB_s2_3|...|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|allele1A_allele1B|0|0|0|0|0|0|...|0|0|0|0|0|0|...|
|allele2A_allele2B|0|0|0|0|0|0|...|0|0|0|0|0|0|...|
|allele3A_allele3B|780.324|906.261|796.327|0|0|0|...|712.676|832.739|733.673|0|0|0|...|
|...|||||||||||||||

SA_s1_1: The first repeat (1) of the allele of subgenome A (SA) in sample1 (s1).

SB_s2_3: The third repeat (3) of the allele of subgenome B (SB) in sample2 (s2).

**Here is one example (The count values of allele pairs in the three replicates of the first sample.):**

```
cd ./00_data/RNA_seq
cat allelepairs.count_selected.txt
# Zo_SA_chr01A_1-Zo_SB_chr01B_47	0	0	0	0	0	0
# Zo_SA_chr01A_2-Zo_SB_chr01B_48	828	907	757	354	384	406
# Zo_SA_chr01A_3-Zo_SB_chr01B_49	0	0	0	0	0	0
# Zo_SA_chr01A_4-Zo_SB_chr01B_50	2619	2758	2141	3305	4140	3122
# Zo_SA_chr01A_5-Zo_SB_chr01B_51	15	21	10	0	0	0
```

- **Output files**
**The expected output files are written in different folders. Here, we have** “pairs.genelist”**,** “name_list.txt”**,** “AvsB.”**,** “stats_number/”**,** **and** “High_expression/” **in the** “03_DEG/1_Class_alleles” **folder, and we have** “.name”**,** “.tpm”**,** “tpm.box.txt”**,** “tpm_boxplot.pdf” **in the** “03_DEG/2_Diff_comparison” **folder, and also** “all.kaks” **and** “.pdf”**, in the “KaKs” folder.**

```
cd 03_DEG/1_Class_alleles
```

|Output files|Description|
|---|---|
|pairs.genelist|Allele pairs|
|name_list.txt|Name of samples|
|AvsB.|Differentially expressed alleles in all samples of subgenomes A and B, respectively (“up” represents alleles that are highly expressed in subgenome A compared to subgenome B.)|
|stats_number/|Multiple statistics of allele expression in different tissues or treatments|
|High_expression/|The number of highly expressed alleles in each chromosome|

```
cd 03_DEG/2_Diff_comparison
cd TPM
```

|Output files|Description|
|---|---|
|.name|The pairs of differentially expressed alleles in Diff0, Diff2 and Diff8|
|.tpm|TPM values of differentially expressed allele pairs in Diff0, Diff2 and Diff8|
|tpm.box.txt|Summary of tpm values for the differentially expressed alleles in each group (Diff0, Diff2 and Diff8)|
|stats_number/|Multiple statistics of allele expression in different tissues or treatments|
|tpm_boxplot.pdf|Distribution diagram of alleles with different differential expression folds under each group (Diff0, Diff2 and Diff8)|

```
cd KaKs
```

|Output files|Description|
|---|---|
|all.kaks|Summary of Ka, Ks and Ka/Ks values between allele pairs in each group (Diff0, Diff2 and Diff8)|
|.pdf|Distribution diagram of Ka, Ks and Ka/Ks between allele pairs in  each group (Diff0, Diff2 and Diff8)|

References to different file formats in bioinformatics:

fasta: https://en.wikipedia.org/wiki/FASTA_format 

gff: https://en.wikipedia.org/wiki/General_feature_format 

bed: https://en.wikipedia.org/wiki/BED_(file_format)

## 5 Parameter setting

### 5.1 Parameter setting in the allele identification step

```
bash ./bin/allele_identification.sh -p chrpairs.txt -a SA -b SB
```

|Parameter|Description|
|---|---|
|-p|Necessary parameter. the target chromosome list|
|-a|Necessary parameter. name of subgenome A|
|-b|Necessary parameter. name of subgenome B.|
|-h|Print brief help message|

**Note:** The names of the subgenomes and the prefixes of the input files must be the same. (eg: If -a is SA, the input file will be SA.pep; if -a is SubgenomeA, the input file will be SubgenomeA.pep.)

The parameter setting requirements in 5.2 are the same.

### 5.2 Parameter setting in the allele specific expression step

```
bash ./bin/allele_specific_expression.sh -a SA -b SB -c allelepairs.count_selected.txt -t allelepairs.tpm_selected.txt -s 102
```

|Parameter|Description|
|---|---|
|-a|Necessary parameter. name of subgenome A|
|-b|Necessary parameter. name of subgenome B|
|-c|Necessary parameter. count matrix of allele pairs|
|-t|Necessary parameter. tpm matrix of allele pairs|
|-s|Necessary parameter. number of samples|
|-h|Print brief help message|

## 6 Note

In the subcommands tpm_boxplot.r (line 27) and kaks_boxplot.r (lines 28,43 and 58) of allele_specific_expression.sh, users can adjust according to their actual results if necessary.


## Reference
Lehmann, R. (2013). 3σ-rule for outlier detection from the viewpoint of geodetic adjustment. Journal of Surveying Engineering, 139(4), 157-165. doi:10.1061/(ASCE)SU.1943-5428.0000112

Sun, P., et al. (2022). WGDI: A user-friendly toolkit for evolutionary analyses of whole-genome duplications and ancestral karyotypes. Molecular Plant, 15(12), 1841-1851. doi:10.1016/j.molp.2022.10.018

Chen, Y., et al. (2020). A collinearity-incorporating homology inference strategy for connecting emerging assemblies in the triticeae tribe as a pilot practice in the plant pangenomic era. Molecular Plant, 13(12), 1694-1708. doi:10.1016/j.molp.2020.09.019







