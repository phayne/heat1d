#!/bin/bash
# script to get ltim_avg values for Diviner mapping cycles
# 9 Oct 2019

#listfile="div_mapcycles_RA.txt"
listfile="cycles2.txt"
while IFS= read cycle; do
  year=${cycle::4}
  STfile='dgdr_st_clc_cyl_'${cycle}'n_128_img.img'
  lbfile='dgdr_st_clc_cyl_'${cycle}'n_128_img.lbl'
  lblink='https://pds-geosciences.wustl.edu/lro/lro-l-dlre-4-rdr-v1/lrodlr_1001/data/gdr_l3/'${year}'/cylindrical/img/'${lbfile}
  STlink='https://pds-geosciences.wustl.edu/lro/lro-l-dlre-4-rdr-v1/lrodlr_1001/data/gdr_l3/'${year}'/cylindrical/img/'${STfile}
  wget -O ${STfile} ${STlink}
  wget -O ${lbfile} ${lblink}
done < "$listfile"
