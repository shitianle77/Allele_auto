#!/bin/bash

######################################################
#Parameter Passing
while getopts "a:b:c:t:s:h" opt; do
  case $opt in
    a)
      SA=$OPTARG
      ;;
    b)
      SB=$OPTARG
      ;;
    c)
      count=$OPTARG
      ;;
    t)
      tpm=$OPTARG
      ;;
    s)
      sample_num=$OPTARG
      ;;
    h)
     echo ""
      ;;
#    \?)
#     echo "Invalid option: -$OPTARG"
#     ;;
  esac
done
######################################################
#Help documentation

 display_usage() {
         echo -e "\n\tThe script is used to study the differential expression of alleles\n"
         echo -e "\tUsage:  bash allele_specific_expression.sh -a [] -b [] -c [] -t [] -s []"
                 echo -e "\tExample:bash allele_specific_expression.sh -a SA -b SB -c allelepairs.count_selected.txt -t allelepairs.tpm_selected.txt -s 102\n"
         echo -e "\t-a: Necessary parameter. name of subgenome A\n"
         echo -e "\t-b: Necessary parameter. name of subgenome B\n"
         echo -e "\t-c: Necessary parameter. count matrix of allele pairs\n"
         echo -e "\t-t: Necessary parameter. tpm matrix of allele pairs\n"
         echo -e "\t-s: Necessary parameter. number of samples\n"
                 echo -e "\t-h or --help: Usage\n"
         }


 # if less than one arguments supplied, display usage
         if [  $# -le 0 ]
         then
                 display_usage
                 exit 1
         fi

 # check whether user had supplied -h or --help . If yes display usage
         if [[ ( $* == "--help") ||  $* == "-h" ]]
         then
                 display_usage
                 exit 0
         fi

######################################################
#Variable Determination

 if [ ! $SA ]|| [ ! $SB ]|| [ ! $count ]|| [ ! $tpm ]|| [ ! $sample_num ];then
         echo -e "\n\tERROR: Missing necessary files."
         exit 1
 fi

#######################################################
#Work Path
pwd=`pwd`

#######################################################
#PartⅡ Allele specific expression
#######################################################
echo -e "\n\t##########################################\n\t##   PartⅡ Allele specific expression   ##\n\t##########################################\n" && sleep 2s

#######################################################
#1、Data processing

cd $pwd/00_data/RNA_seq

########################################################
#2、Identification and classification of differentially expressed allele pairs
########################################################
echo -e "\n\t##  1. Identification and classification of differentially expressed allele pairs  ##\n" && sleep 2s

########################################################
cd $pwd
mkdir 03_DEG
cd 03_DEG

mkdir 1_Class_alleles
cd 1_Class_alleles

##run
Rscript $pwd/bin/allele_DEseq.r ./ $pwd/00_data/RNA_seq/${count} ${sample_num}

##stats
echo a*|sed 's/ /\n/g'|sed 's/\./\t/'|cut -f 1|sort|uniq |sed 's/a/a\t/g'|sort -k2,2n|sed 's/a\t/a/g' > name_list.txt
cut -f 1 $pwd/00_data/RNA_seq/${count} > pairs.genelist

for i in `cat name_list.txt`
do 
cat ${i}.diff00.txt ${i}.diff2.txt ${i}.diff0.txt ${i}.diff8.txt|grep -v 'baseMean'|cut -f 1 > ${i}_right.genelist
cat pairs.genelist ${i}_right.genelist |sort|uniq -c|awk '{print $1}'|sort|uniq -c > ${i}_stats
cat pairs.genelist ${i}_right.genelist |sort|uniq -c|awk '{if($1==1)print $2}' > ${i}_NA.txt
done

########################################################
#3、Comparative analysis between subgenomes
########################################################
echo -e "\n\t##               2. Comparative analysis between subgenomes               ##\n" && sleep 2s

#################################################
#3.1 Differences in the expression of alleles in different tissues or treatments 
#################################################
echo -e "\n\t# Differences in the expression of alleles in different tissues or treatments #\n" && sleep 2s

#################################################
#Result integration
mkdir stats_number
mv *_stats ./stats_number
cd stats_number

##Noexpression
for i in `cat ../name_list.txt`
do
awk -va=${i} '{print a"\t"$1"\t"$2}' ${i}_stats|awk -va="Noexpression" '{if($3=="1")print $1"\t"$2"\t"a}'
done > stats_Noexpression
##Diff00
for i in `cat ../name_list.txt`
do
grep -v 'baseMean' ../${i}.diff00.txt |wc -l |awk -va=${i} -vb="diff00" '{print a"\t"$1"\t"b}' 
done > stats_diff00
##Diff0
for i in `cat ../name_list.txt`
do
grep -v 'baseMean' ../${i}.diff0.txt |wc -l |awk -va=${i} -vb="diff0" '{print a"\t"$1"\t"b}' 
done > stats_diff0
##Diff2
for i in `cat ../name_list.txt`
do
grep -v 'baseMean' ../${i}.diff2.txt |wc -l |awk -va=${i} -vb="diff2" '{print a"\t"$1"\t"b}' 
done > stats_diff2
##Diff8
for i in `cat ../name_list.txt`
do
grep -v 'baseMean' ../${i}.diff8.txt |wc -l |awk -va=${i} -vb="diff8" '{print a"\t"$1"\t"b}' 
done > stats_diff8

cat stats_Noexpression stats_diff00 stats_diff0 stats_diff2 stats_diff8 |sed 's/a/a\t/'|sort -k2.2n |sed 's/a\t/a/'| awk '{print $2"\t"$1"\t"$3}' > duidietu.txt
rm a* stats_Noexpression stats_diff00 stats_diff0 stats_diff2 stats_diff8

# Plot
cd $pwd/03_DEG/1_Class_alleles/stats_number

Rscript $pwd/bin/Classification_histogram.r ./ duidietu.txt ../name_list.txt

#################################################
# 3.2 The number of highly expressed alleles in each chromosome
#################################################
echo -e "\n\t# The number of highly expressed alleles in each chromosome #\n" && sleep 2s

#################################################
#3.2.1 Data processing
cd $pwd/03_DEG/1_Class_alleles
mkdir High_expression
cd High_expression

for i in `cat ../name_list.txt`
do
cat ../${i}.diff0.sup.txt |awk 'NR>1{print $0}'|sed 's/-/\t/'|awk '{print $1}'|grep -f - ../../../00_data/${SA}.gff|cut -f 1|sort |uniq -c|awk '{print $1"\t"$2}'|sed 's/A/\tA/'|awk '{print $1"\t"$3"\t"$2}'
done >> high_expression_A.txt

for i in `cat ../name_list.txt`
do
cat ../${i}.diff0.dom.txt |awk 'NR>1{print $0}'|sed 's/-/\t/'|awk '{print $2}'|grep -f - ../../../00_data/${SB}.gff|cut -f 1|sort |uniq -c|awk '{print $1"\t"$2}'|sed 's/B/\tB/'|awk '{print $1"\t"$3"\t"$2}'
done >> high_expression_B.txt

cat high_expression_A.txt high_expression_B.txt|sed '1iNumber\tGroup\tChromosomes' > Highexpression_chr_number.txt
rm high_expression_A.txt high_expression_B.txt

#3.2.2 Plot
cd $pwd/03_DEG/1_Class_alleles/High_expression

Rscript $pwd/bin/Highexpression_chr_number.r ./ Highexpression_chr_number.txt

#################################################
# 3.3 The distribution of alleles with different differential expression folds under each condition
#################################################
echo -e "\n\t# The distribution of alleles with different differential expression folds under each condition #\n" && sleep 2s

#################################################
cd $pwd/03_DEG/
mkdir 2_Diff_comparison
cd 2_Diff_comparison

##（1）tpm
mkdir TPM
cd TPM

cat $pwd/03_DEG/1_Class_alleles/*.diff0.txt|grep -v 'baseMean'|cut -f 1|sort|uniq > diff0.name
cat $pwd/03_DEG/1_Class_alleles/*.diff2.txt|grep -v 'baseMean'|cut -f 1|sort|uniq > diff2.name
cat $pwd/03_DEG/1_Class_alleles/*.diff8.txt|grep -v 'baseMean'|cut -f 1|sort|uniq > diff8.name
grep -f diff0.name $pwd/00_data/RNA_seq/allelepairs.tpm_selected.txt > diff0.tpm
grep -f diff2.name $pwd/00_data/RNA_seq/allelepairs.tpm_selected.txt > diff2.tpm
grep -f diff8.name $pwd/00_data/RNA_seq/allelepairs.tpm_selected.txt > diff8.tpm

cut -f 2- diff0.tpm  | sed 's/\t/\n/g' | awk '{print $1 "\tdiff0"}' >  tpm.box.tpm
cut -f 2- diff2.tpm  | sed 's/\t/\n/g' | awk '{print $1 "\tdiff2"}' >> tpm.box.tpm
cut -f 2- diff8.tpm  | sed 's/\t/\n/g' | awk '{print $1 "\tdiff8"}' >> tpm.box.tpm
grep -v 'e' tpm.box.tpm > tpm.box.txt

##Plot
cd $pwd/03_DEG/2_Diff_comparison/TPM

Rscript $pwd/bin/tpm_boxplot.r ./ tpm.box.txt

##（2）ka
cd $pwd/03_DEG/2_Diff_comparison
mkdir KaKs
cd KaKs

awk '{print $1"-"$2}' $pwd/02_wgdi/allele_pairs_lasted.txt > allelepairs_lasted.txt
awk 'NR>1{print $1"-"$2"\t"$3"\t"$4}' $pwd/02_wgdi/SA_SB.ks.txt |awk '{if($3==0)print $1"\t"$2"\t"$3"\t"0; else print $1"\t"$2"\t"$3"\t"($2/$3)}' |grep -f allelepairs_lasted.txt - > kaks.txt

grep -f $pwd/03_DEG/2_Diff_comparison/TPM/diff0.name kaks.txt|awk -va="diff0" '{print $0"\t"a}' > diff0.kaks
grep -f $pwd/03_DEG/2_Diff_comparison/TPM/diff2.name kaks.txt|awk -va="diff2" '{print $0"\t"a}' > diff2.kaks
grep -f $pwd/03_DEG/2_Diff_comparison/TPM/diff8.name kaks.txt|awk -va="diff8" '{print $0"\t"a}' > diff8.kaks

cat  diff0.kaks diff2.kaks diff8.kaks|grep -v 'e' > all.kaks
rm kaks.txt diff0.kaks diff2.kaks diff8.kaks

##Plot
cd $pwd/03_DEG/2_Diff_comparison/KaKs
Rscript $pwd/bin/kaks_boxplot.r ./ all.kaks

#rm $pwd/03_DEG/1_Class_alleles/a*
