#!/bin/sh -f

# run python scripts: runPTCLM.py to configure/build/run a case
# 
# NOTE: (1) user must provide non-default input datasets as indiciate in the following
#       (2) this DEMO configure case in /lustre/pfs1/cades-ccsi/proj-shared/project_acme, and,
#           run case in /lustre/pfs1/cades-ccsi/scratch/$USER

#0.1 case setup directory: this script where hold, and could be directory to hold 'cases'
CASESETUP_DIR=${PWD}
CCSI_DIR=/lustre/pfs1/cades-ccsi

CLM_PF_TOOLS=$CCSI_DIR/proj-shared/models/alm-pf-tools

#1.1 codes and default input data direcotries
CLM_COUPLED_MODEL=$CCSI_DIR/proj-shared/models/ACME-fmyuan

#1.2 user-defined root working directory (scratch) 
WORKDIR=$CCSI_DIR/scratch/$USER/
mkdir -p ${WORKDIR}

MYINPUTDATA=$CCSI_DIR/proj-shared/project_acme/ACME_inputdata

#2.1 create case/run directory, if not yet
CLM_CASE_DIR=$CCSI_DIR/proj-shared/project_acme/cases
mkdir -p ${CLM_CASE_DIR}
CLM_RUN_DIR=${WORKDIR}/runs
mkdir -p ${CLM_RUN_DIR}

#2.2 Run python scripts to create/build/submit a run. Please refer to the python scripts 'runPTCLM.py' with --help
cd ${CLM_PF_TOOLS}
python ./runPTCLM.py --site=US-Brw \
                 --sitegroup=AmeriFlux \
                 --caseidprefix=alm \
                 --caseroot=${CLM_CASE_DIR} \
                 --runroot=${CLM_RUN_DIR} \
                 --ccsm_input=${MYINPUTDATA} \
                 --cesmdir=${CLM_COUPLED_MODEL} \
                 --compset=I1850CLM45CN \
                 --coldstart \
                 --ad_spinup \
                 --CH4 \
                 --nitrif_denitrif \
                 --no_fire \
                 --restart_n=1 --restart_units=nyears \
                 --run_n=10 --run_units=nyears \
                 --machine=cades --mach_dir=${CLM_PF_TOOLS}/machines-userdefined \
                 --compiler=gnu \
                 --mpilib=openmpi \
                 --rmold \
                 --clean_config \
                 --clean_build \
                 --debug_build \
                 --xpts=1 \
                 --ypts=1 \
                 --update-datm-domain \
                 --hist_userdefined=clm_output_adspinup \
                 --np=1 --ppn=1 \
                 --surfdatafile=surfdata_1x1pt_US-Brw_simyr1850.nc \
                 --walltime=01:00:00 \
                 --jobname=testalm

cd ${CASESETUP_DIR}

##----------------------------#
## this is for point/multiple-points(regional) ALM and alm-pflotran model setup/runs on CADES-CCSI
## Please direct your comments, questions or ask for instructions to
## 
## Fengming YUAN, CCSI/ESD-ORNL
## (865)57408486(O)
## yuanf@ornl.gov
## 
## 2016-05-25
##----------------------------#
