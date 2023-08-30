#!/bin/bash

######################################################
#Parameter Passing
while getopts "p:i:h" opt; do
  case $opt in
    p)
      chrpairs=$OPTARG
      ;;
    i)
      IQR=$OPTARG
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
         echo -e "\n\tThis script is for fiiltering allele pairs\n"
         echo -e "\tUsage:  bash allele_fiiltering.sh -p [] -i []"
                 echo -e "\tExample:bash allele_fiiltering.sh -p chrpairs.txt -i 1.5\n"
         echo -e "\t-p: Necessary parameter. homologous chromosome pairs\n"
         echo -e "\t-i: Necessary parameter. the multiples of the IQR (interquartile range)\n"
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

 if [ ! $chrpairs ]|| [ ! $IQR ];then
         echo -e "\n\tERROR: Missing necessary files."
         exit 1
 fi
#######################################################
#Work Path
pwd=`pwd`

#######################################################
#3.3 Filtering
echo -e "\n\t# Filtering-Tukey #\n" && sleep 2s

#######################################################
cd $pwd
mkdir 02_wgdi-Tukey && cd 02_wgdi-Tukey

#data processing
cp $pwd/02_wgdi/genepairs_info.tsv ./
cp $pwd/02_wgdi/block.tsv ./

cat genepairs_info.tsv|awk -F"\t" '{OFS="\t"}NR>1{x=($5-$4)/2+$4;y=($9-$8)/2+$8;slope=(y/x); print $0"\t"x"\t"y"\t"slope}' |sed '1iBlock\tgeneID1\tchr1\tstart_pos1\tend_pos1\tgeneID2\tchr2\tstart_pos2\tend_pos2\tka_NG86\tks_NG86\tka_YN00\tks_YN00\tvalue_x\tvalue_y\tslope' > slope_genepairs_info.tsv

#run(filtering)
Rscript $pwd/bin/Filtering_Tukey.r ./ genepairs_info.tsv ${IQR} ks_NG86
Rscript $pwd/bin/Filtering_Tukey.r ./ slope_genepairs_info.tsv ${IQR} slope

#
cat ks_NG86-${IQR}.genepairs_info.tsv|awk 'NR>1{if($14!="NA")print $0}'|cut -f 1-13 > aa
grep 'Block' ks_NG86-${IQR}.genepairs_info.tsv |cat - aa|cut -f 1-13 > filtered_ks_NG86-${IQR}.genepairs_info.tsv
rm aa

cat slope-${IQR}.genepairs_info.tsv|awk 'NR>1{if($17!="NA")print $0}'|cut -f 1-13 > aa
grep 'Block' slope-${IQR}.genepairs_info.tsv |cat - aa|cut -f 1-13 > filtered_slope-${IQR}.genepairs_info.tsv
rm aa slope_genepairs_info.tsv

cat filtered_ks_NG86-${IQR}.genepairs_info.tsv filtered_slope-${IQR}.genepairs_info.tsv|grep -v 'Block'|sort|uniq -c|awk '{if($1=="2")print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}'|sort -k1,1n|sed '1iBlock\tgeneID1\tchr1\tstart_pos1\tend_pos1\tgeneID2\tchr2\tstart_pos2\tend_pos2\tka_NG86\tks_NG86\tka_YN00\tks_YN00' > raw_filtered_ks_NG86.genepairs_info.tsv 

#Obtain filtered Co-linear block and allele pair files
awk 'NR>1{print $1}' filtered_ks_NG86-${IQR}.genepairs_info.tsv|sort|uniq|sort -k1,1n > keep.ksblock.id
awk 'NR>1{print $1}' filtered_slope-${IQR}.genepairs_info.tsv|sort|uniq|sort -k1,1n > keep.positionblock.id

cat keep.ksblock.id keep.positionblock.id|sort|uniq -c |awk '{if($1==2){print $2}}'|sed '1i Block' |sort -k1,1 |awk 'NR==FNR{a[$1]=$0}NR>FNR{if($1 in a)print $0}' - raw_filtered_ks_NG86.genepairs_info.tsv > filtered.ks-slope-bioplot_genepairs_info.tsv
rm keep.ksblock.id keep.positionblock.id

#Final allele pairs obtained
cat filtered.ks-slope-bioplot_genepairs_info.tsv|awk 'NR>1{print $2"\t"$6}' > allele_pairs_lasted.txt

#Collinearity of a pair of alleles on two subgenomes
cd $pwd/02_wgdi-Tukey
mkdir allele_plot && cd allele_plot
# run
for i in `cat $pwd/00_data/$chrpairs|cut -f 3`
do
awk 'NR==FNR{a[$4]=$0;}NR!=FNR{print $0,a[$1]}' $pwd/01_genetribe/${i}/${i}A.bed $pwd/02_wgdi-Tukey/allele_pairs_lasted.txt|awk '{if(NF=="8")print $0}' > ${i}.aa
awk 'NR==FNR{a[$4]=$0;}NR!=FNR{print $0,a[$2]}' $pwd/01_genetribe/${i}/${i}B.bed ${i}.aa|awk -va=${i} '{print a"\t"$4"\t"$5"\t"$10"\t"$11}' > ${i}.coord.allele
done

rm *.aa raw_filtered_ks_NG86.genepairs_info.tsv ks_NG86-1.5.genepairs_info.tsv slope-1.5.genepairs_info.tsv
cat *.coord.allele > pairs.coord.allele

#plot
Rscript $pwd/bin/plot_pairs_allele.r ./ pairs.coord.allele

