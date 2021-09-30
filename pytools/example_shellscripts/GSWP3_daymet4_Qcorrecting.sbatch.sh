#!/bin/sh

#  GSWP3_daymet4_Qcorrecting.sbatch.sh
#  
#
#  Created by Yuan, Fengming on 10/12/21.
#  

#SBATCH --time=01:00:00
#SBATCH -A ccsi
#SBATCH -p batch
#SBATCH --mem=125G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node 1
#SBATCH --job-name=merge_tiledata
#SBATCH -o ./%j-output.txt
#SBATCH -e ./%j-error.txt
#SBATCH  --exclusive

source $MODULESHOME/init/bash

#module unload python
#module load python/3.9.4

cwd=$(pwd)

cd ${cwd}

GSWP3_daymet=/nfs/data/ccsi/f9y/GSWP3_daymet
TILE_NO=TILE14236

E3SM_INPUTDATA=/lustre/or-scratch/cades-ccsi/proj-shared/project_acme/e3sm_inputdata

python3 ./metdata_daymet_Qcorrecting.py \
--mettype='cplbypass_GSWP3_daymet4' \
--e3sm_metdomain=$GSWP3_daymet/$TILE_NO/cpl_bypass_full/domain.nc \
--e3sm_metdir=$GSWP3_daymet/$TILE_NO/cpl_bypass_full \
--user_mettype='cplbypass_GSWP3' \
--user_metdir=$E3SM_INPUTDATA/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v2.c180716/cpl_bypass_full \
--varname='QBOT'  --nc_write --ncout_mettype='cplbypass_GSWP3_daymet4' \
--ncout_stdmetdir=$GSWP3_daymet/$TILE_NO/cpl_bypass_full

echo "done!"

cd ${cwd}
