

#1. ######prepare gff file
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_all_16_modify/state_variability/
ls *.bed | while read id;
do
    echo $id
    cat $id | cut -f 1-2 | sed 's/\t/:/g' > 1.txt
    cat $id | cut -f 3 | sed 's/^/-/g'  > 2.txt
    paste 1.txt 2.txt | sed 's/\t//g' > 3.txt
    paste $id 3.txt | awk '{print $1,$2,$3,$5,$4}' | sed 's/ /\t/g' > AA_state_id/$(basename $id ".bed")_id.bed  
    paste $id 3.txt | awk '{print $1, $5, $5, $2, $3, $3-$2, $7=".",$8="-",$5}' | sed 's/ /\t/g' > AA_super_enhancer/$(basename $id ".bed").gff
    cat $id > AA_super_enhancer/$id

done


cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_all_16_modify/state_variability/AA_super_enhancer
for tissue in Adipose BMarrow Bursa Cecum Cerebellum Colon Cortex Duodenum Gizzard Heart Hypothalamus Ileum Jejunum Kidney Liver Lung Muscle Provent ShellGland Spleen Testis Thymus Trachea
do
   if [ -f ${tissue}_E6.gff ] ; then
   cat ${tissue}_E6.gff ${tissue}_E7.gff ${tissue}_E8.gff ${tissue}_E9.gff ${tissue}_E10.gff > ${tissue}_enhancer.gff

   fi
done




#2. ######prepare bam file
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_all_16_modify/state_variability/AA_super_enhancer

ls /group/zhougrp/zhangyuan/Chicken/Aligned_Reads/*_CA.bam | while read id; do (ln -s $id $(basename $id)); done
ls /group/zhougrp/zhangyuan/Chicken/Aligned_Reads/*_CB.bam | while read id; do (ln -s $id $(basename $id)); done
ls /group/zhougrp/zhangyuan/Chicken/Aligned_Reads/*bam.bai| while read id; do (ln -s $id $(basename $id)); done






#3. creat refseq.ucsc

#!/bin/bash
SPECIES=$1   ##hg38
FEATURE=$2   ##gene transcipt
GTFFILE=$3   ##/group/zhougrp/zhangyuan/genome/human/GRCh38_Ensembl100.gtf
mkdir -p annotation
echo -e "#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\tX\tX\tX\t\tX\tname2" > annotation/$SPECIES"_refseq.ucsc"

if [[ $FEATURE == "gene" ]]; then
awk -F'[\t ]' '{
  if($3=="gene")
    print "0\t" $10 "\tchr" $1 "\t" $7 "\t" $4 "\t" $5 "\t" $4 "\t" $5 "\t.\t.\t.\t.\t" $10}' $GTFFILE | sed s/\"//g | sed 's/;//g' >> annotation/$SPECIES"_refseq.ucsc"

elif [[ $FEATURE == "transcript" ]]; then
awk -F'[\t ]' '{
  if($3=="transcript")
    print "0\t" $10 "\tchr" $1 "\t" $7 "\t" $4 "\t" $5 "\t" $4 "\t" $5 "\t.\t.\t.\t.\t" $10}' $GTFFILE | sed s/\"//g | sed 's/;//g' >> annotation/$SPECIES"_refseq.ucsc"
fi
echo "Annotation file: "$SPECIES"_refseq.ucsc"


cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_all_16_modify/state_variability/AA_super_enhancer
sbatch -p high -t 1-0 -c 1 ~/software/ROSE/refseq.sh galGal6 gene /group/zhougrp/zhangyuan/genome/chicken/galGal6_Ensembl100.gtf






#4. modify genomeDict
vim /home/zhypan/software/ROSE/bin/ROSE_main.py
    genomeDict = {
        'HG18':'%s/annotation/hg18_refseq.ucsc' % (cwd),
        'MM9': '%s/annotation/mm9_refseq.ucsc' % (cwd),
        'HG19':'%s/annotation/hg19_refseq.ucsc' % (cwd),
        'HG38':'%s/annotation/hg38_refseq.ucsc' % (cwd),
        'MM8': '%s/annotation/mm8_refseq.ucsc' % (cwd),
        'MM10':'%s/annotation/mm10_refseq.ucsc' % (cwd),
        'GALGAL6':'%s/annotation/galGal6_refseq.ucsc' % (cwd),
        }





#5. identify the super enhancer

#!/bin/bash
PATHTO=/home/zhypan/software/ROSE
PYTHONPATH=$PATHTO/lib
export PYTHONPATH
export PATH=$PATH:$PATHTO/bin
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_all_16_modify/state_variability/AA_super_enhancer
for tissue in Adipose BMarrow Bursa Cecum Cerebellum Colon Cortex Duodenum Gizzard Heart Hypothalamus Ileum Jejunum Kidney Liver Lung Muscle Provent ShellGland Spleen Testis Thymus Trachea
do
echo $tissue
ROSE_main.py -g galGal6 -i ${tissue}_enhancer.gff -r H3K27ac_${tissue}_CA.bam -c Input_${tissue}_CA.bam -o AA_super_result/${tissue}_enhancer_CA -s 12500 -t 2500
ROSE_main.py -g galGal6 -i ${tissue}_enhancer.gff -r H3K27ac_${tissue}_CB.bam -c Input_${tissue}_CB.bam -o AA_super_result/${tissue}_enhancer_CB -s 12500 -t 2500

done


sbatch -p bmh -t 2-0 -c 24 --mem 48G 114_super_enhaner_chicken.sh 












