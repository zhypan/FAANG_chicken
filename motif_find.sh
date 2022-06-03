
#!/bin/sh
PATH=$PATH:/home/ckern/homer/bin

ls /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_all_16_modify/state_variability/AARegulatory_module/AA_TSR_E6/*.bed | while read id;
do
out=/group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_all_16_modify/state_variability/AARegulatory_module/AA_TSR_E6/AA_Motif_anohter/$(basename $id ".bed")_motif
echo $out
findMotifsGenome.pl $id /group/zhougrp/zhangyuan/genome/chicken/galGal6.fa $out -len 8,10,12 -size 200 -mask -p 24
done