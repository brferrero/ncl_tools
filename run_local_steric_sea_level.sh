#!/usr/bin/env bash
#
# execute local_steric_sea_level.ncl script
#  

VARIN="RHO"
#SCENARIO="B20TRC5CNBDRD"
SCENARIO="BRCP85C5CNBDRD"
REFYEAR=2006

DIRIN="/chuva/wk2/bruno/cesm_lens_slr/${SCENARIO}/"${VARIN}"/"
FILELIST=$(ls ${DIRIN})
VAROUT="ssl"
DIROUT="./nc/"

for FIN in $FILELIST; do
    echo "File: "$FIN
    
    FILEIN=${DIRIN}${FIN}
    FOUT=$(echo ${FIN} | sed 's/'${VARIN}'/'${VAROUT}'/g')
    FILEOUT=${DIROUT}${FOUT}
    ncl filein=\"${FILEIN}\" fileout=\"${FILEOUT}\" varin=\"${VARIN}\" varout=\"${VAROUT}\" refyear=${REFYEAR} local_steric_sea_level.ncl
done

