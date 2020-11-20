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
  --varname=TBOT --nc_create --ncout_standard --ncout_cplbypass

echo "PRECT"
python3 ./metdata_merge.py --mettype='Site,ATS_h5' \
  --e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015 \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --user_metdir=/Users/f9y/ATS_ROOT/INPUT_DATA \
  --user_metfile='CESM-RCP8.5-2006-2100_dm1985-2015.h5'  \
  --varname='PRECT' --source_adjusting --nc_write --ncout_standard --ncout_cplbypass

echo "PSRF"
python3 ./metdata_merge.py --mettype='Site,ATS_h5' \
  --e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015 \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --user_metdir=/Users/f9y/ATS_ROOT/INPUT_DATA \
  --user_metfile='CESM-RCP8.5-2006-2100_dm1985-2015.h5'  \
  --varname='PSRF' --nc_write --ncout_standard --ncout_cplbypass

echo "RH"
python3 ./metdata_merge.py --mettype='Site,ATS_h5' \
  --e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015 \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --user_metdir=/Users/f9y/ATS_ROOT/INPUT_DATA \
  --user_metfile='CESM-RCP8.5-2006-2100_dm1985-2015.h5'  \
  --varname='RH' --nc_write --ncout_standard --ncout_cplbypass

echo "WIND"
python3 ./metdata_merge.py --mettype='Site,ATS_h5' \
  --e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015 \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --user_metdir=/Users/f9y/ATS_ROOT/INPUT_DATA \
  --user_metfile='CESM-RCP8.5-2006-2100_dm1985-2015.h5'  \
  --varname='WIND' --nc_write --ncout_standard --ncout_cplbypass

echo "FSDS"
python3 ./metdata_merge.py --mettype='Site,ATS_h5' \
  --e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015 \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --user_metdir=/Users/f9y/ATS_ROOT/INPUT_DATA \
  --user_metfile='CESM-RCP8.5-2006-2100_dm1985-2015.h5'  \
  --varname='FSDS' --nc_write --ncout_standard --ncout_cplbypass

echo "FLDS"
python3 ./metdata_merge.py --mettype='Site,ATS_h5' \
  --e3sm_metdir=/Users/f9y/clm4_5_inputdata/atm/datm7/CLM1PT_data/1x1pt_US-Brw_19852015 \
  --e3sm_metheader='????-??' \
  --e3sm_metdomain=/Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_Grid.nc \
  --user_metdir=/Users/f9y/ATS_ROOT/INPUT_DATA \
  --user_metfile='CESM-RCP8.5-2006-2100_dm1985-2015.h5'  \
  --varname='estFLDS' --nc_write --ncout_standard --ncout_cplbypass

echo "DONE!"




