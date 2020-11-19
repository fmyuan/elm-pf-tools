#!/bin/bash

source ~/.bashrc

# data sub-daily scaling and merging
echo "TBOT"
python3 ./metdata_merge.py --mettype='Site,ATS_h5' \
  --e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015 \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --user_metdir=/Users/f9y/ATS_ROOT/INPUT_DATA \
  --user_metfile='CESM-RCP8.5-2006-2100_dm1985-2015.h5'  \
  --varname=TBOT --nc_create

echo "PRECT"
python3 ./metdata_merge.py --mettype='Site,ATS_h5' \
  --e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015 \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --user_metdir=/Users/f9y/ATS_ROOT/INPUT_DATA \
  --user_metfile='CESM-RCP8.5-2006-2100_dm1985-2015.h5'  \
  --varname='PRECT' --source_adjusting --nc_write

echo "PSRF"
python3 ./metdata_merge.py --mettype='Site,ATS_h5' \
  --e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015 \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --user_metdir=/Users/f9y/ATS_ROOT/INPUT_DATA \
  --user_metfile='CESM-RCP8.5-2006-2100_dm1985-2015.h5'  \
  --varname='PSRF' --nc_write

echo "RH"
python3 ./metdata_merge.py --mettype='Site,ATS_h5' \
  --e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015 \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --user_metdir=/Users/f9y/ATS_ROOT/INPUT_DATA \
  --user_metfile='CESM-RCP8.5-2006-2100_dm1985-2015.h5'  \
  --varname='RH' --nc_write

echo "WIND"
python3 ./metdata_merge.py --mettype='Site,ATS_h5' \
  --e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015 \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --user_metdir=/Users/f9y/ATS_ROOT/INPUT_DATA \
  --user_metfile='CESM-RCP8.5-2006-2100_dm1985-2015.h5'  \
  --varname='WIND' --nc_write

echo "FSDS"
python3 ./metdata_merge.py --mettype='Site,ATS_h5' \
  --e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015 \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --user_metdir=/Users/f9y/ATS_ROOT/INPUT_DATA \
  --user_metfile='CESM-RCP8.5-2006-2100_dm1985-2015.h5'  \
  --varname='FSDS' --nc_write

echo "FLDS"
python3 ./metdata_merge.py --mettype='Site,ATS_h5' \
  --e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015 \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --user_metdir=/Users/f9y/ATS_ROOT/INPUT_DATA \
  --user_metfile='CESM-RCP8.5-2006-2100_dm1985-2015.h5'  \
  --varname='estFLDS' --nc_write

echo "DONE!"




