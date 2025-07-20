#!/bin/bash

set -e

# Created 2025-03-10 10:05:06

CASEDIR="/gpfs/wolf2/cades/cli185/proj-shared/$USER/project-e3sm/cases/ELMngee4_pan-arctic-CAVM0.5deg_GSWP3_ICB1850CNRDCTCBC_ad_spinup"

/gpfs/wolf2/cades/cli185/proj-shared/f9y/models/E3SM_ORNL_IM/cime/scripts/create_newcase --case "${CASEDIR}" --mach cades-baseline --compset ICB1850CNRDCTCBC --res ELM_USRDAT --mpilib openmpi-amanzitpls --walltime 24:0:00 --handle-preexisting-dirs u --project e3sm --compiler gnu

cd "${CASEDIR}"

./xmlchange MOSART_MODE=NULL

./xmlchange DIN_LOC_ROOT=/gpfs/wolf2/cades/cli185/world-shared/e3sm/inputdata

./xmlchange DIN_LOC_ROOT_CLMFORC=/gpfs/wolf2/cades/cli185/world-shared/e3sm/inputdata/atm/datm7/

./xmlchange ELM_USRDAT_NAME=pan-arctic_CAVM.0.5deg

./xmlchange --append ELM_BLDNML_OPTS="-bgc_spinup on"

./xmlchange ATM_DOMAIN_PATH=/gpfs/wolf2/cades/cli185/world-shared/e3sm/inputdata/share/domains/domain.clm

./xmlchange LND_DOMAIN_PATH=/gpfs/wolf2/cades/cli185/world-shared/e3sm/inputdata/share/domains/domain.clm

./xmlchange ATM_DOMAIN_FILE=domain.lnd.r05_RRSwISC6to18E3r5.240328_cavm1d.nc

./xmlchange LND_DOMAIN_FILE=domain.lnd.r05_RRSwISC6to18E3r5.240328_cavm1d.nc

./xmlchange DOUT_S=FALSE

./xmlchange ATM_NCPL=24

./xmlchange NTASKS=1280

./xmlchange NTHRDS=1

./xmlchange MAX_TASKS_PER_NODE=128

./xmlchange MAX_MPITASKS_PER_NODE=128

./xmlchange STOP_OPTION=nyears

./xmlchange STOP_N=200

./xmlchange REST_N=20

./xmlchange --id PIO_TYPENAME --val netcdf


echo "
&clm_inparm
 hist_mfilt = 1
 hist_nhtfrq = 0
 hist_empty_htapes = .true.
 hist_fincl1 = 'TLAI','TOTSOMC','TSOI_10CM','SOILWATER_10CM'
 hist_dov2xy = .true.
 fsurdat = '/gpfs/wolf2/cades/cli185/world-shared/e3sm/inputdata/lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr1850_c240308_TOP_cavm1d.nc'
 stream_fldfilename_ndep = '/gpfs/wolf2/cades/cli185/world-shared/e3sm/inputdata/lnd/clm2/ndepdata/fndep_elm_cbgc_exp_simyr1849-2101_1.9x2.5_ssp245_c240903.nc'
 nyears_ad_carbon_only = 25
 spinup_mortality_factor = 10
 metdata_type = 'gswp3'
 metdata_bypass = '/gpfs/wolf2/cades/cli185/world-shared/e3sm/inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v2.c180716/cpl_bypass_full'
 co2_file = '/gpfs/wolf2/cades/cli185/world-shared/e3sm/inputdata/atm/datm7/CO2/fco2_datm_rcp4.5_1765-2500_c130312.nc'
 aero_file = '/gpfs/wolf2/cades/cli185/world-shared/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_rcp4.5_monthly_1849-2104_1.9x2.5_c100402.nc'
">>user_nl_elm

./case.setup

echo "
string(APPEND CPPDEFS " -DCPL_BYPASS")
">>cmake_macros/universal.cmake


./case.build

#./case.submit


