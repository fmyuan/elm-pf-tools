#!/bin/bash -f
#SBATCH --time=2:00:00
#SBATCH -A ccsi
#SBATCH -p burst
#SBATCH --mem=0G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node 1
#SBATCH --job-name=utar_daymet
#SBATCH -o ./%j-output.txt
#SBATCH -e ./%j-error.txt
#SBATCH  --exclusive 

source $MODULESHOME/init/bash

cwd=$(pwd)

#####
DAYMET_SRCFILE=/lustre/or-scratch/cades-ccsi/scratch/f9y/Daymet_ESM 
TARNO=149
####XX=28
for XX in {90..93}

do
  TILE=${TARNO}${XX}

  if tar -xf ${DAYMET_SRCFILE}_${TARNO}XX.tar --wildcards --no-anchored "*/${TILE}/*.nc" &
  sleep 10
  then
    wait

    echo "DONE untaring with TILE ${TILE} !"

  else
    echo "ERROR untaring with TILE ${TILE} !"
    exit &?
  fi

  ####
  bash ./mv_years_daymetfiles.sh ${TILE}

  echo "done!"

done

cd ${cwd}

