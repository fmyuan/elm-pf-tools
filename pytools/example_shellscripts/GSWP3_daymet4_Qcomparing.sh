#!/bin/bash -f

source ~/.bashrc
cwd=$(pwd)

python3 ../metdata_daymet_Qcorrecting.py \
--e3sm_metdomain=/Users/f9y/clm4_5_inputdata/atm/datm7/GSWP3_daymet/TILE14236/cpl_bypass_full_ToolikFS/domain.nc \
--mettype='cplbypass_GSWP3_daymet4' \
--e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/GSWP3_daymet/TILE14236/cpl_bypass_full \
--user_mettype='cplbypass_GSWP3' \
--user_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v2.c180716_ToolikFS-Grid/cpl_bypass_full \
--varname='QBOT'  --nc_create \
--ncout_mettype='cplbypass_GSWP3_daymet4' \
--ncout_stdmetdir=/Users/f9y/clm4_5_inputdata/atm/datm7/GSWP3_daymet/TILE14236/cpl_bypass_full/corrected
