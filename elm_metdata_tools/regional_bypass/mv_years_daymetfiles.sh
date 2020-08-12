#!/bin/bash -f

cwd=$(pwd)

DAYMET_ROOT=/nfs/data/ccsi/f9y/GSWP3_daymet

TILE=$1 ###e.g. 13867 #daymet tile index

cd ./

#move data into one directory
mkdir -p ${DAYMET_ROOT}/TILE${TILE}
mkdir -p ${DAYMET_ROOT}/TILE${TILE}/Precip3Hrly
mkdir -p ${DAYMET_ROOT}/TILE${TILE}/Solar3Hrly
mkdir -p ${DAYMET_ROOT}/TILE${TILE}/TPHWL3Hrly
mkdir -p ${DAYMET_ROOT}/TILE${TILE}/cpl_bypass_full

for IY in {1980..2014}
do
  if mv "./${IY}/${TILE}/Precip3Hrly"/*.nc "${DAYMET_ROOT}/TILE${TILE}/Precip3Hrly/"
  then
    wait  
  else
    exit $?
  fi
  
  if mv "./${IY}/${TILE}/Solar3Hrly"/*.nc "${DAYMET_ROOT}/TILE${TILE}/Solar3Hrly/"
  then
    wait
  else
    exit $?
  fi
  
  if mv "./${IY}/${TILE}/TPHWL3Hrly"/*.nc "${DAYMET_ROOT}/TILE${TILE}/TPHWL3Hrly/"
  then
    wait
  else
    exit $?
  fi
  
  rm -rf ./${IY}/
done

wait
echo "DONE!"


