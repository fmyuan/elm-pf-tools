#!/bin/bash

source ~/.bashrc

# checking after sub-daily scaling and merging
echo "TBOT"
python3 ./metdata_merge.py --mettype='Site' \
  --e3sm_metdir=./ \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --varname=TBOT --plotting

echo "PRECT"
python3 ./metdata_merge.py --mettype='Site' \
  --e3sm_metdir=./ \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --varname='PRECT' --plotting

echo "PSRF"
python3 ./metdata_merge.py --mettype='Site' \
  --e3sm_metdir=./ \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --varname='PSRF' --plotting

echo "RH"
python3 ./metdata_merge.py --mettype='Site' \
  --e3sm_metdir=./ \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --varname='RH' --plotting

echo "WIND"
python3 ./metdata_merge.py --mettype='Site' \
  --e3sm_metdir=./ \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --varname='WIND' --plotting

echo "FSDS"
python3 ./metdata_merge.py --mettype='Site' \
  --e3sm_metdir=./ \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --varname='FSDS' --plotting

echo "FLDS"
python3 ./metdata_merge.py --mettype='Site' \
  --e3sm_metdir=./ \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --varname='FLDS' --plotting

echo "DONE!"




