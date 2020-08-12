#!/bin/bash -f
#SBATCH --time=24:00:00
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
DAYMET_SRCFILE=Daymet_ESM_138XX.tar
### For AK SP: DAYMET tile index - 13867~13870
TILE=13867

cd ${cwd}

### extract daymet forcing data from GAO's *.tar files


if tar -xf ${DAYMET_ROOT}/${DAYMET_SRCFILE} --wildcards --no-anchored "*/${TILE}/*.nc" &
sleep 10
then
  wait

  echo "DONE untaring with TILE ${TILE} !"

else
  exit &?
fi

cd ${cwd}

