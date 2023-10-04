#!/bin/bash -f

source ~/.bashrc
cwd=$(pwd)

python3 ../metdata_comparing.py \
--metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_ToolikFS-GRID_navy.nc \
--mettype='Site' \
--metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_ToolikFS --metheader='????-??' \
--mettype2='cplbypass_GSWP3' \
--metdir2=/Users/f9y/clm4_5_inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v2.c180716_ToolikFS-Grid/cpl_bypass_full \
--user_mettype='cplbypass_GSWP3_daymet' \
--user_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/GSWP3_daymet/TILE14236 \
--varname='QBOT'  --plotting --similarity --yrmin=1989 --yrmax=2014

