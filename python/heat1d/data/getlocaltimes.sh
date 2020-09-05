#!/bin/bash
# script to get ltim_avg values for Diviner mapping cycles
# 9 Oct 2019

latmin="20.0"
latmax="21.0"
lonmin="-49.0"
lonmax="-48.0"

ppd="128"

#declare -a cycles=("20090705" "20090713" "20090810")
#for i in "${cycles[@]}"

infile="div_mapcycles_RA.txt"
echo $infile
while IFS= read cycle; do
  echo $cycle
  outfile="LT_${cycle}_bin.out"
  echo $outfile
  divgdr debug=2 res=${ppd} date=${cycle} dn=N extract=ltim_avg \
  | pcons lat=${latmin},${latmax} lon=${lonmin},${lonmax} > ${outfile}
done < "$infile"
