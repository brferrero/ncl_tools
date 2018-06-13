#!/usr/bin/env bash
#
# execute volume_wgt_mean.ncl script
#  
# BASIN={all full atl ind pac}

VARIN="TEMP"
DIRIN="/chuva/wk2/bruno/cesm_lens_slr/B20TRC5CNBDRD/"${VARIN}"/"
FIN="b.e11.B20TRC5CNBDRD.f09_g16.008.pop.h."${VARIN}".yr.192001-200512.nc"
FILEIN=${DIRIN}${FIN}

BASIN="all"
VAROUT="gmtemp"
DIROUT="./"
FOUT="b.e11.B20TRC5CNBDRD.f09_g16.008.pop.h."${VAROUT}"."${BASIN}".yr.192001-200512.nc"
FILEOUT=${DIROUT}${FOUT}

ncl filein=\"${FILEIN}\" fileout=\"${FILEOUT}\" varin=\"${VARIN}\" varout=\"${VAROUT}\" basin=\"${BASIN}\" volume_wgt_mean.ncl
