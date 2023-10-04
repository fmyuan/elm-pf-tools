#!/bin/bash -f

source ~/.bashrc
cwd=$(pwd)

python3 ../Daymet_tilemerge.py \
  --daymet_elm_mapfile=daymet_elm_mappings.txt \
  --workdir='../TILE11208/cpl_bypass_full-KXTN' \
  --workdir2='../TILE11209/cpl_bypass_full-KXTN,../TILE11388/cpl_bypass_full-KXTN,../TILE11389/cpl_bypass_full-KXTN' \
  --mapfile_redoxy
