#!/bin/sh -f

# run python scripts: runPTCLM.py to configure/build/run a case
# 
# NOTE: (1) user must provide non-default gcc/mpich/netcdf as indiciate in the 'env_mach_specific.mymac'
#       (2) this DEMO configure case in $HOME/project_acme/cases, and,
#           run case in $HOME/project_acme/scratch/runs, and,
#       (3) inputdata in $HOME/project_acme/cesm-inputdata

#0.1 case setup directory: this script where hold, and could be directory to hold 'cases'
CASESETUP_DIR=${PWD}

CLM_PF_TOOLS=$HOME/project_acme/alm-pf-tools

#1.1 codes
CLM_COUPLED_MODEL=$HOME/mygithub/ACME-fmyuan

#1.2 user-defined root working directory (scratch) 
WORKDIR=$HOME/project_acme/scratch
mkdir -p ${WORKDIR}

#1.3 inputdata
MYINPUTDATA=$HOME/project_acme/cesm-inputdata

#2.1 create case/run directory, if not yet
CLM_CASE_DIR=$HOME/project_acme/cases
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
                 --machine=mymac --mach_dir=${CLM_PF_TOOLS}/machines-userdefined \
                 --compiler=gnu \
                 --mpilib=mpich \
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
## this is for point/multiple-points(regional) ALM and alm-pflotran model setup/runs on Mac OS
## with user-built GCC/MPICH and libraries such as hdf5/netcdf/petsc specified in 'env_mach_specific.mymac'
## Please direct your comments, questions or ask for detail instructions to
## 
## Fengming YUAN, CCSI/ESD-ORNL
## (865)57408486(O)
## yuanf@ornl.gov
## 
## 2016-05-25
##----------------------------#
