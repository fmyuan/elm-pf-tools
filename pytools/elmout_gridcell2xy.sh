#!/bin/bash -f
#SBATCH --time=2:00:00
#SBATCH -A cli185
#SBATCH -p batch
#SBATCH --mem=0G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node 1
#SBATCH --job-name=gridcell2xy
#SBATCH -o ./%j-output.txt
#SBATCH -e ./%j-error.txt
#SBATCH  --exclusive 

source ~/.bashrc

cwd=$(pwd)

cd ${cwd}

srun -n 1 python3 /gpfs/wolf2/cades/cli185/proj-shared/f9y/models/elm-pf-tools/pytools/ELMout_daymet_tile.py \
--elmheader=/gpfs/wolf2/cades/cli185/proj-shared/f9y/archives/elm_SPAK2/ELM2021_AKSP_1km_ICB20TRCNPRDCTCBC0/run/ELM2021_AKSP_1km_ICB20TRCNPRDCTCBC0.elm.h1.199 \
--elm_varname='ALL' \
--daymet_elm_mapfile=projection_elm_mappings.txt \
--proj_name='AlaskaAber' --mapfile_redoxy

