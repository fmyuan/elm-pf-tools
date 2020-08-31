#!/bin/bash -f

cwd=$(pwd)

DAYMET_ROOT=/nfs/data/ccsi/f9y/GSWP3_daymet

VARNAME=$1 ###e.g. 13867 #variable name
echo "${VARNAME}"

cd ./

YRS_ZONE=1980-2014_z01
if [ -f "./temp_all.nc" ]; then
  rm ./temp_all.nc
fi

if [ -f "./zone_mappings.txt" ]; then
  rm ./zone_mappings.txt
fi

for TILE in {14047..14050}
do
  NCFILE=${DAYMET_ROOT}/TILE${TILE}/cpl_bypass_full/GSWP3_daymet4_${VARNAME}_${YRS_ZONE}.nc
  # make gridcell dimension 'n' as record dimension so that can be concatated
  if ncks -h -O --mk_rec_dmn n ${NCFILE} temp${TILE}.nc
  then
    wait  
  else
    exit $?
  fi

  #concatate nc file by gridcell dimension 'n'
  if [ -f "./temp_all.nc" ]; then 
    ncrcat -O -h ./temp_all.nc ./temp${TILE}.nc -o ./temp_all.nc
    
  else
    cp ./temp${TILE}.nc ./temp_all.nc
  fi

  #
  if cat ${DAYMET_ROOT}/TILE${TILE}/cpl_bypass_full/zone_mappings.txt >>./zone_mappings.txt
  then
    wait
  else
    exit $?
  fi
    
  #rm -rf ./${NCFILE}/
done

#fix gridcell dimension 'n'
if ncks -h -O --fix_rec_dmn n temp_all.nc -o temp_all.nc
then
  wait  
else
  exit $?
fi


wait
echo "DONE!"


