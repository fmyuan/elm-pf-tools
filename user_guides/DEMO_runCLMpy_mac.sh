#!/bin/sh

# run python scripts: runCLM.py to build/run a case
# on MAC : eg. $bash DEMO_runCLMpy_mac.sh
# NOTE: user must provide non-default input datasets as indiciate in the following

#0.1 case setup directory: this script where hold, ususally up directory of /clm4-pf-tools
CASESETUP_DIR=${PWD}
CLM4_PF_TOOLS_DIR = CASESETUP_DIR/clm4-pf-tools
#CLM4_PF_TOOLS_DIR =/projects/cesm/bitbucket_copy/clm85-pflotran-ngee-sci

#1.1 codes and default input data direcotries - REQUIRE editing by USER
CLM_COUPLED_MODEL=/projects/cesm/bitbucket_copy/clm85-pflotran-ngee-sci
CESMINPUT0=/projects/cesm/clm45inputdata0
USERINPUT=/projects/cesm/cesminputdata

#1.2 change to user-defined working directory for inputdata/cases/runs - REQUIRE editng by USER
WORKDIR=${HOME}/ngee_simulations

MYINPUTDATA=${WORKDIR}/myinputdata
$${CLM4_PF_TOOLS_DIR}/link_dirtree ${CESMINPUT0} ${MYINPUTDATA}

#2.1 DATM domain and data used actually is for 1x1pt
cd ${MYINPUTDATA}/atm/datm7/domain.clm
rm domain.lnd.1x1pt_US-Brw_navy.nc
ln -sf ${USERINPUT}/atm/datm7/domain.clm/domain.lnd.1x1pt_US-Brw_navy.nc domain.lnd.1x1pt_US-Brw_navy.nc

cd ${MYINPUTDATA}/atm/datm7
rm -rf 1x1pt_US-Brw
ln -sf ${USERINPUT}/atm/datm7/1x1pt_US-Brw 1x1pt_US-Brw

#2.2 Need to use user-defined CLM surfdata and its domain
cd ${MYINPUTDATA}/lnd/clm2/surfdata_map
rm -f surfdata_1x1pt_US-Brw_simyr1850.nc
rm -r surfdata_1x1pt_US-Brw_simyr2000.nc
ln -sf ${USERINPUT}/lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850_ngeeC_default.nc ./surfdata_1x1pt_US-Brw_simyr1850.nc
ln -sf ${USERINPUT}/lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850_ngeeC_default.nc ./surfdata_1x1pt_US-Brw_simyr2000.nc

cd ${MYINPUTDATA}/share/domains/domain.clm
ln -sf ${USERINPUT}/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy_ngeeC.nc domain.lnd.1x1pt_US-Brw_navy.nc

#3.1 create case/run directory, if not yet
CLM_CASE_DIR=${WORKDIR}/cases
mkdir -p ${CLM_CASE_DIR}
CLM_RUN_DIR=${WORKDIR}/runs
mkdir -p ${CLM_RUN_DIR}

#3.2 Run python scripts to create/build/submit a run. Please refer to the python scripts 'runCLM.py' with --help
cd ${CLM4_PF_TOOLS_DIR}
python ./runCLM.py --site=US-Brw \
                 --sitegroup=AmeriFlux \
                 --caseidprefix=clmr85_1x1pt \
                 --caseroot=${CLM_CASE_DIR} \
                 --runroot=${CLM_RUN_DIR} \
                 --ccsm_input=${MYINPUTDATA} \
                 --cesmdir=${CLM_COUPLED_MODEL} \
                 --compset=I1850CLM45CN \
                 --coldstart \
                 --ad_spinup \
                 --vertsoilc \
                 --CH4 \
                 --nitrif_denitrif \
                 --no_fire \
                 --nyears_ad_spinup=10 \
                 --machine=userdefined \
                 --osname=Darwin \
                 --compiler=gnu \
                 --mpilib=mpich \
                 --rmold \
                 --clean_config \
                 --clean_build \
                 --debug_build \
                 --nopointdata \
                 --xpts=1 \
                 --ypts=1 \
                 --update-datm-domain \
                 --hist_userdefined=clm-output-adspinup \
                 --np=1 \
                 --ninst=1 \
                 --jobname=run_stdclm --no_submit

cd ${CASESETUP_DIR}

