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
### For AK SP: DAYMET tile index - 13867~13870
###TILE=13870

cd ${cwd}

### For AK SP: DAYMET tile index - 13867~13870

### For AK above N60: 

#DAYMET_SRCFILE=Daymet_ESM_135XX.tar
#for TILE in {13507..13520}

#DAYMET_SRCFILE=Daymet_ESM_136XX.tar
#for TILE in {13685..13699}

#DAYMET_SRCFILE=Daymet_ESM_137XX.tar
#for TILE in {13700..13700}

DAYMET_SRCFILE=Daymet_ESM_138XX.tar
for TILE in {13866..13880}

#DAYMET_SRCFILE=Daymet_ESM_140XX.tar
#for TILE in {14047..14060}

#DAYMET_SRCFILE=Daymet_ESM_142XX.tar
#for TILE in {14227..14240}

#DAYMET_SRCFILE=Daymet_ESM_144XX.tar
#for TILE in {14409..14419}


do

  ### extract 1 tile of daymet forcing data from GAO's *.tar files


  if tar -xf ${DAYMET_ROOT}/${DAYMET_SRCFILE} --wildcards --no-anchored "*/${TILE}/*.nc" &
  sleep 10
  then
    wait

    echo "DONE untaring with TILE ${TILE} !"

  else
    echo "ERROR untaring with TILE ${TILE} !"
    exit &?
  fi

done


cd ${cwd}

