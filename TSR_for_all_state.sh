#!/bin/bash
####get each state for each tissues
mkdir state_variability
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_16_modify
ls *16_segments.bed | while read id;
do
  echo $id
  b=$(basename $id "_16_segments.bed")
  echo $b
  for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
  do
     grep -w $state $id | sort -k1,1 -k2,2n > "state_variability/"$b"_"$state".bed"

  done
done  


#Gs
####get Gs: total regions of each state (Gs) across 14 tissues 
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_16_modify/state_variability
mkdir AAGs
for i in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
cat *$i".bed" |sort -k1,1 -k2,2n > 1.bed
bedtools merge -i 1.bed > "AAGs/"$i"_Gs.bed"
done

#!/bin/bash
mkdir AARegulatory_module
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_16_modify/state_variability
ls *${state}*bed | while read id;
do
echo $id
bedtools intersect -a <(sort -k1,1 -k2,2n /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_16_modify/state_variability/AAGs/${state}_Gs.bed) -b <(sort -k1,1 -k2,2n ${id}) -c -sorted >  "AARegulatory_module/"${id%%.*}"_Gs.bed"
done


echo chr start end > 1.txt
ls *_${state}.bed >> 1.txt
cat 1.txt | cut -d "_" -f 1 | perl -p -e 's/\n/ /g'| sed '$ s/.$/\n/' > AARegulatory_module/header.txt
#paste all tissues together
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_16_modify/state_variability/AARegulatory_module
ls *${state}_Gs.bed  | while read id;
do
  echo $id
  cat ${id} | cut -f 4 | paste -s >> 2.txt
done
awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}' 2.txt > 3.txt
rm 2.txt
cut  -f 1-3 Cecum_${state}_Gs.bed > 5.txt
paste 5.txt 3.txt |sed 's/ /\t/g' >6.txt
cat header.txt 6.txt |sed 's/ /\t/g' > all_${state}_Gs.csv




#normalized the one count and count the number for each region
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_16_modify/state_variability/AARegulatory_module
cat all_${state}_Gs.csv | cut -f 4-18 |  sed 's/20/1/g' | sed 's/21/1/g'| sed 's/22/1/g'| sed 's/23/1/g'|sed 's/24/1/g'|sed 's/25/1/g'|sed 's/10/1/g'| sed 's/11/1/g'| sed 's/12/1/g'| sed 's/13/1/g'| sed 's/14/1/g'| sed 's/15/1/g' | sed 's/16/1/g'| sed 's/17/1/g'|  sed 's/18/1/g'| sed 's/19/1/g'| sed 's/2/1/g'| sed 's/3/1/g'| sed 's/4/1/g'| sed 's/5/1/g' | sed 's/6/1/g'| sed 's/7/1/g'|  sed 's/8/1/g'| sed 's/9/1/g' > 3.txt
cat 3.txt | awk '{for(i=1;i<=NF;i++){a[NR]+=$i}print $0,a[NR]}' > 4.txt
cut  -f 1-3 all_${state}_Gs.csv > 5.txt
paste  5.txt 4.txt |sed 's/ /\t/g' > all_${state}_Gs_one_count.csv



for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
echo $state
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_all_16_modify/state_variability/AARegulatory_module
mkdir AA_TSR_${state}
cat all_${state}_Gs_one_count.csv | awk '{if (($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==23) print $0}' > AA_TSR_${state}/TSR_All_common.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($7+$9+$11+$15+$16)==5 && ($4+$5+$6+$8+$10+$12+$13+$14+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Intestine_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($8+$10+$14)==3 && ($4+$5+$6+$7+$9+$11+$12+$13+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Brain_common_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if (($5+$6+$23+$26)==4 && ($4+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$24+$25)==0) print $0}' > AA_TSR_${state}/TSR_Immune_common_${state}.txt

cat all_${state}_Gs_one_count.csv | awk '{if ($4==1 && ($5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Adipose_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($5==1 && ($4+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_BMarrow_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($6==1 && ($4+$5+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Bursa_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($7==1 && ($4+$5+$6+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Cecum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($8==1 && ($4+$5+$6+$7+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Cerebellum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($9==1 && ($4+$5+$6+$7+$8+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Colon_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($10==1 && ($4+$5+$6+$7+$8+$9+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Cortex_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($11==1 && ($4+$5+$6+$7+$8+$9+$10+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Duodenum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($12==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Gizzard_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($13==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Heart_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($14==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Hypothalamus_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($15==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Ileum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($16==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Jejunum_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($17==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$18+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Kidney_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($18==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$19+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Liver_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($19==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$20+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Lung_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($20==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$21+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Muscle_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($21==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$22+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Provent_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($22==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$23+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_ShellGland_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($23==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Spleen_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($24==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Testis_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($25==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$26)==0) print $0}' > AA_TSR_${state}/TSR_Thymus_${state}.txt
cat all_${state}_Gs_one_count.csv | awk '{if ($26==1 && ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25)==0) print $0}' > AA_TSR_${state}/TSR_Trachea_${state}.txt

done


for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
cat all_${state}_Gs_one_count.csv | awk '{if (($6+$23)==3 && ($4+$5+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$24+$25+$26)==0) print $0}' > AA_TSR_${state}/TSR_Immune_common_${state}.txt

done


cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_all_16_modify/state_variability/AARegulatory_module
for state in E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14 E15
do
rm AA_TSR_${state}/TSR_Immune_common_${state}.txt
rm AA_TSR_${state}/TSR_Immune_common_${state}_id.bed
done