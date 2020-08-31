#!/bin/bash -f
#SBATCH --time=48:00:00
#SBATCH -A ccsi
#SBATCH -p batch
#SBATCH --mem=125G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node 1
#SBATCH --job-name=untar_daymet
#SBATCH -o ./%j-output.txt
#SBATCH -e ./%j-error.txt
#SBATCH  --exclusive 

source $MODULESHOME/init/bash

cwd=$(pwd)

DAYMET_ROOT=/nfs/data/ccsi/f9y/GSWP3_daymet

cd ${cwd}

### For AK above N60: 

DAYMET_SRCFILE=Daymet_ESM_138XX.tar
$YEAR=1980
tar -cf AK-${DAYMET_SRCFILE} ./1980/

for YEAR in {1981..2014}
do
  if tar --append --file=AK-${DAYMET_SRCFILE} ./${YEAR}/ &
  sleep 10
  then
    wait
    echo "DONE tar  ${YEAR} with AK-${DAYMET_SRCFILE} !"
  else
    echo "ERROR tar  ${YEAR} with AK-${DAYMET_SRCFILE} !"
    exit &?
  fi

done

cd ${cwd}

