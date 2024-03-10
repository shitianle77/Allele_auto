#!/bin/bash

######################################################
#Parameter Passing
while getopts "p:a:b:h" opt; do
  case $opt in
    p)
      chrpairs=$OPTARG
      ;;
    a)
      SA=$OPTARG
      ;;
    b)
      SB=$OPTARG
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
         echo -e "\n\tThis script is for identifing allele pairs\n"
         echo -e "\tUsage:  bash allele_identification.sh -p [] -a [] -b []"
                 echo -e "\tExample:bash allele_identification.sh -p chrpairs.txt -a SA -b SB\n"
         echo -e "\t-c: Necessary parameter. the target chromosome list\n"
         echo -e "\t-a: Necessary parameter. name of subgenome A\n"
         echo -e "\t-b: Necessary parameter. name of subgenome B\n"
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

 if [ ! $chrpairs ]|| [ ! $SA ]|| [ ! $SB ];then
         echo -e "\n\tERROR: Missing necessary files."
         exit 1
 fi
#######################################################
#Work Path
pwd=`pwd`

#######################################################
#PartⅠ Allele identification
#######################################################
echo -e "\n\t##########################################\n\t##     Part Ⅰ Allele identification     ##\n\t##########################################\n" && sleep 2s

#######################################################
#1、Data processing
cd $pwd/
mkdir 01_genetribe
cd 01_genetribe

#（1）bed information for each chromosome
for i in `cat $pwd/00_data/$chrpairs|cut -f 3` 
do 
mkdir ${i} && cd ${i}
grep ${i} $pwd/00_data/${SA}.bed > ${i}A.bed
grep ${i} $pwd/00_data/${SB}.bed > ${i}B.bed
cd ..
done

#（2）protein sequences for each chromosome
for i in `cat $pwd/00_data/$chrpairs|cut -f 3|sed 's/chr//'`
do
cd chr${i}
perl $pwd/bin/fasta_no_blank.pl $pwd/00_data/${SA}.pep|grep -A 1 ${i}A > chr${i}A.fa
perl $pwd/bin/fasta_no_blank.pl $pwd/00_data/${SB}.pep|grep -A 1 ${i}B > chr${i}B.fa
cd ..
done

#（3）chromosome group information
for i in `cat $pwd/00_data/$chrpairs|cut -f 3`
do
cd ${i}
echo 'chrNA' > ${i}A.chrlist
echo 'chrNB' > ${i}B.chrlist
cd ..
done

########################################################
#2、Identification of homologous genes between homologous chromosomes
## genetribe-run
########################################################
echo -e "\n\t##  1. Identification of homologous genes between homologous chromosomes  ##\n" && sleep 2s

#######################################################
cd $pwd/01_genetribe

#（1）Preparing parallelisation commands
for i in `cat $pwd/00_data/$chrpairs|cut -f 3`
do
echo "cd ${i}; genetribe core -l ${i}A -f ${i}B"
done > genetribe.txt

#（2）Running (parallelisation)
cd $pwd/01_genetribe
ParaFly -c genetribe.txt -CPU 8

#（3）Data integration
cat ./*/*.RBH > raw_RBH.genepairs
cat ./*/*.RBH |awk '{print $1"\t"$2}' > RBH.genepairs

########################################################
#3、Filteri
## WGDI-dependent
########################################################
echo -e "\n\t##                              2. Filtering                              ##\n" && sleep 2s

########################################################
cd $pwd/01_genetribe
mkdir stats
cd stats

#######################################################
#3.1 Extraction of blast results for homologous gene pairs
echo -e "\n\t# Extraction of blast results for homologous gene pairs #\n" && sleep 2s

#######################################################
for i in `cat $pwd/00_data/$chrpairs|cut -f 3`
do
ln -s ../${i}/output/${i}A_${i}B.blast2 ./
ln -s ../${i}/output/${i}B_${i}A.blast2 ./
ln -s ../${i}/${i}A_${i}B.RBH ./
awk '{print $2"\t"$1}' ${i}B_${i}A.blast2 |cat - ${i}A_${i}B.blast2|awk '{print $1"\t"$2}'|sort|uniq -c|awk '{if($1==2)print $2"\t"$3}'|grep -w -f - ../RBH.genepairs > ${i}_RBH_2hits.genepairs2 
grep -w -f ${i}_RBH_2hits.genepairs2 ${i}A_${i}B.blast2 > ${i}_RBH_2hits.blast2
awk '{print $1"\t"$2}' ${i}A_${i}B.RBH|cat - ${i}_RBH_2hits.genepairs2|sort|uniq -c|awk '{if($1==1)print $2"\t"$3}'|grep -w -f - ${i}A_${i}B.blast2 > ${i}_RBH_1hits_A2B.blast2
awk '{print $1"\t"$2}' ${i}A_${i}B.RBH|cat - ${i}_RBH_2hits.genepairs2|sort|uniq -c|awk '{if($1==1)print $2"\t"$3}'|grep -w -f - ${i}A_${i}B.blast2|awk '{print $1"\t"$2}' > ${i}_RBH_1hits_A2B.genepairs2
awk '{print $1"\t"$2}' ${i}A_${i}B.RBH|cat - ${i}_RBH_2hits.genepairs2 ${i}_RBH_1hits_A2B.genepairs2|sort|uniq -c|awk '{if($1==1)print $3"\t"$2}'|grep -w -f - ${i}B_${i}A.blast2|awk '{print $2"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}'|cat ${i}_RBH_2hits.blast2 ${i}_RBH_1hits_A2B.blast2 - > ${i}_RBH.blast1
rm *.RBH *.blast2 *.genepairs2
done

#Data integration
cat *.blast1 > ../SA_SB.blast
cd ..
rm -rf stats

for i in `cat $pwd/00_data/$chrpairs|cut -f 3`
do
cd ${i}
rm -rf output
cd ..
done

#######################################################
#3.2 Calculate Co-linear blocks and ks
echo -e "\n\t# Calculate Co-linear blocks and ks #\n" && sleep 2s

#######################################################
cd $pwd
mkdir 02_wgdi
cd 02_wgdi

#3.2.1 Data processing
#(1) gff
ln -s $pwd/01_genetribe/SA_SB.blast ./

awk '{print $1"\t"$2}' SA_SB.blast |sed 's/\t/\n/' > SA_SB.genelist
grep -w -f SA_SB.genelist $pwd/00_data/${SA}.gff > SA.gff
grep -w -f SA_SB.genelist $pwd/00_data/${SB}.gff > SB.gff

## pep and cds sequences
cat $pwd/00_data/${SA}.cds| sed 'N; s/\n/ /'|grep -w -f SA_SB.genelist - |sed 's/ /\n/' > SA.cds.fa
cat $pwd/00_data/${SB}.cds| sed 'N; s/\n/ /'|grep -w -f SA_SB.genelist - |sed 's/ /\n/' > SB.cds.fa
cat $pwd/00_data/${SA}.pep| sed 'N; s/\n/ /'|grep -w -f SA_SB.genelist - |sed 's/ /\n/' > SA.pep.fa
cat $pwd/00_data/${SB}.pep| sed 'N; s/\n/ /'|grep -w -f SA_SB.genelist - |sed 's/ /\n/' > SB.pep.fa

cat SA.cds.fa SB.cds.fa > SA_SB.cds.fa
cat SA.pep.fa SB.pep.fa > SA_SB.pep.fa
rm  SA_SB.genelist SA.cds.fa SB.cds.fa SA.pep.fa SB.pep.fa

#3.2.2 run
cd $pwd/02_wgdi
source /home/ytbao/bin/miniconda3/bin/activate wgdi

ln -s $pwd/bin/conf/SA_SB_collinearity.conf ./
ln -s $pwd/00_data/${SA}.len ./SA.len
ln -s $pwd/00_data/${SB}.len ./SB.len
wgdi -icl SA_SB_collinearity.conf

ln -s $pwd/bin/conf/SA_SB_ks.conf ./
wgdi -ks SA_SB_ks.conf

ln -s $pwd/bin/conf/SA_SB_block.conf ./
wgdi -bi SA_SB_block.conf

#ln -s $pwd/bin/conf/SA_SB_blockks.conf ./
#sed -e "4s/SA/${SA}/" -e "5s/SB/${SB}/" SA_SB_collinearity.conf -i
#wgdi -bk SA_SB_blockks.conf

#######################################################
#3.3 Filtering
echo -e "\n\t# Filtering #\n" && sleep 2s

#######################################################
cd $pwd/02_wgdi

#data processing
cat SA_SB_block.csv|sed 's/,/\t/g' > block.tsv
awk '{print NF"\t"$0}' SA_SB.ks.txt|awk '{if($1=="2")print $0"\t"0"\t"0"\t"0"\t"0; else print $0}'|cut -f 2-|awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($2 in a)print a[$2]"\t"$0}' - SA.gff |awk '{print $2"\t"$1"\t"$7"\t"$9"\t"$10"\t"$3"\t"$4"\t"$5"\t"$6}' > aa
awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($2 in a)print a[$2]"\t"$0}' aa SB.gff |awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1"\t"$10"\t"$12"\t"$13"\t"$6"\t"$7"\t"$8"\t"$9}' |sed '1igeneID1\tchr1\tstart_pos1\tend_pos1\tgeneID2\tchr2\tstart_pos2\tend_pos2\tka_NG86\tks_NG86\tka_YN00\tks_YN00' > raw_genepairs_info.tsv
awk '{if($0~/.*#/){split($0,a,"Alignment");print a[1]}else{print $0"\t"a[2]}}'  SA_SB.collinearity.txt |sed '/#/d' |sed 's/ /\t/g'|awk '{print $6"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}'|sed 's/://' > collinearity.txt
awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($2 in a)print a[$2]"\t"$0}' raw_genepairs_info.tsv collinearity.txt |awk '{print $13"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}'|sed '1iBlock\tgeneID1\tchr1\tstart_pos1\tend_pos1\tgeneID2\tchr2\tstart_pos2\tend_pos2\tka_NG86\tks_NG86\tka_YN00\tks_YN00' > genepairs_info.tsv
rm aa raw_genepairs_info.tsv collinearity.txt

pairfile=`ls genepairs_info.tsv`
blockfile=`ls block.tsv`

#3.3.1 Removal of homologous gene pairs with outlying ks values
##Calculate the mean and standard deviation(sd) from the data distribution of pair
mean=`cat ${pairfile}|awk -F "\t" 'NR>1{if($10>0) print $10}'|grep -v "NA"|awk '{x[NR]=$0; s+=$0; n++} END{a=s/n; print a}' `
sd=`cat ${pairfile}|awk -F "\t" 'NR>1{if($10>0) print $10}'|grep -v "NA"|awk '{x[NR]=$0; s+=$0; n++} END{a=s/n; for (i in x){ss += (x[i]-a)^2} sd = sqrt(ss/n); print sd}' `

##Filtering blocks using mean and sd
cat ${blockfile}|awk -vm=${mean} -vs=${sd} -F"\t" '{OFS="\t"}NR>1{plus=m+3*s;minus=m-3*s;if($10>minus && $10<plus) print $1}' > keep.ksblock.id

#3.3.2 Removal of physically located outlier homologous gene pairs
##Calculate the mean and standard deviation(sd) from the data distribution of pair
MEAN=`cat ${pairfile}|awk -F "\t" 'NR>1{x=($5-$4)/2+$4;y=($9-$8)/2+$8;print y/x}'|grep -v "NA"|awk '{x[NR]=$0; s+=$0; n++} END{a=s/n; print a}' `
SD=`cat ${pairfile}|awk -F "\t" 'NR>1{x=($5-$4)/2+$4;y=($9-$8)/2+$8;print y/x}'|grep -v "NA"|awk '{x[NR]=$0; s+=$0; n++} END{a=s/n; for (i in x){ss += (x[i]-a)^2} sd = sqrt(ss/n); print sd}' `

##Filtering blocks using mean and sd
cat ${blockfile}|awk -vm=${MEAN} -vs=${SD} -F"\t" '{OFS="\t"}NR>1{x=($5-$4)/2+$4;y=($7-$6)/2+$6;slope=y/x;plus=m+3*s;minus=m-3*s;if(slope>minus && slope<plus) print $1}' > keep.positionblock.id

#######################################################
#3.4 Obtain filtered Co-linear block and allele pair files
cat keep.ksblock.id keep.positionblock.id|sort|uniq -c |awk '{if($1==2){print $2}}'|sed '1i id' |sort -k1,1 |awk 'NR==FNR{a[$1]=$0}NR>FNR{if($1 in a)print $0}' - ${blockfile} > filtered.${blockfile}
cat keep.ksblock.id keep.positionblock.id|sort|uniq -c |awk '{if($1==2){print $2}}'|sed '1i Block' |sort -k1,1 |awk 'NR==FNR{a[$1]=$0}NR>FNR{if($1 in a)print $0}' - ${pairfile} > filtered.${pairfile}
rm keep.ksblock.id keep.positionblock.id

#######################################################
#4、Final allele pairs obtained
#######################################################
cat filtered.genepairs_info.tsv|awk 'NR>1{print $2"\t"$6}' > allele_pairs_lasted.txt


#######################################################
#5、Collinearity of a pair of alleles on two subgenomes
#######################################################
mkdir allele_plot && cd allele_plot
# run
for i in `cat $pwd/00_data/$chrpairs|cut -f 3`
do
awk 'NR==FNR{a[$4]=$0;}NR!=FNR{print $0,a[$1]}' $pwd/01_genetribe/${i}/${i}A.bed $pwd/02_wgdi/allele_pairs_lasted.txt|awk '{if(NF=="8")print $0}' > ${i}.aa
awk 'NR==FNR{a[$4]=$0;}NR!=FNR{print $0,a[$2]}' $pwd/01_genetribe/${i}/${i}B.bed ${i}.aa|awk -va=${i} '{print a"\t"$4"\t"$5"\t"$10"\t"$11}' > ${i}.coord.allele
done

rm *.aa
cat *.coord.allele > pairs.coord.allele

#plot
Rscript $pwd/bin/plot_pairs_allele.r ./ pairs.coord.allele
