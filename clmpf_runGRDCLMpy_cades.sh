#!/bin/sh -f

# run python scripts: runGRDCLM.py to configure/build/run a global case

# NOTE: (1) user must provide non-default input datasets as indiciate in the following
#       (2) this DEMO configure case in /lustre/pfs1/cades-ccsi/proj-shared/project_acme, and,
#           run case in /lustre/pfs1/cades-ccsi/scratch/$USER

#0.1 case setup directory: this script where hold, and could be directory to hold 'cases'
CASESETUP_DIR=${PWD}

#CCSI_DIR=/lustre/pfs1/cades-ccsi
CCSI_DIR=/lustre/or-hydra/cades-ccsi

#CLM_PF_TOOLS=$CCSI_DIR/proj-shared/models/alm-pf-tools
CLM_PF_TOOLS=./

#1.1 codes and default input data direcotries
CLM_COUPLED_MODEL=$CCSI_DIR/proj-shared/models/ACME-fmyuan

#1.2 user-defined root working directory (scratch) 
WORKDIR=$CCSI_DIR/scratch/$USER/
mkdir -p ${WORKDIR}

MYINPUTDATA=${WORKDIR}/myinputdata
${CLM_COUPLED_MODEL}/scripts/link_dirtree ${CESMINPUT0} ${MYINPUTDATA}

#2.1 DATM domain and data
cd ${MYINPUTDATA}/atm/datm7

# the following is for compset 'I1850CRUCLM45CN', using back-filled CRU/QIAN dataset (1901-1920)
unlink atm_forcing.datm7.cruncep_qianFill.0.5d.V4.c130305
ln -sf ${USERINPUT}/atm/datm7/atm_forcing.datm7.cruncep_qianFill.0.5d.V4.c130305 atm_forcing.datm7.cruncep_qianFill.0.5d.V4.c130305

# the following is for compset 'I1850CLM45CN', using DATM%QIAN (1948 - 1972)
#unlink -rf atm_forcing.datm7.Qian.T62.c080727
#ln -sf ${USERINPUT}/atm/datm7/atm_forcing.datm7.Qian.T62.c080727 atm_forcing.datm7.Qian.T62.c080727

#2.2 Need to use CLM surfdata, maps and its domain
cd ${MYINPUTDATA}/lnd/clm2/surfdata_map
ln -sf ${USERINPUT}/lnd/clm2/surfdata_map/surfdata_1.9x2.5_simyr1850_c130927.nc surfdata_1.9x2.5_simyr1850.nc
ln -sf ${USERINPUT}/lnd/clm2/surfdata_map/surfdata_1.9x2.5_simyr2000_c130927.nc surfdata_1.9x2.5_simyr2000.nc

cd ${MYINPUTDATA}/lnd/clm2/mappingdata/maps
unlink 1.9x2.5
ln -sf ${USERINPUT}/lnd/clm2/mappingdata/maps/1.9x2.5 1.9x2.5

cd ${MYINPUTDATA}/share/domains
ln -sf ${USERINPUT}/share/domains/domain.lnd.fv1.9x2.5_gx1v6.090206_vji.nc domain.lnd.fv1.9x2.5_gx1v6.090206.nc
ln -sf ${USERINPUT}/share/domains/domain.ocn.fv1.9x2.5_gx1v6_090403_vji.nc domain.ocn.fv1.9x2.5_gx1v6_090403.nc

# 2.3 CPL (maps)
cd ${MYINPUTDATA}/cpl/cpl6
unlink map_r05_TO_g16_aave.120920.nc
unlink map_r05_to_gx1v6_e1000r300_090226.nc
ln -sf ${USERINPUT}/cpl/cpl6/map_r05_TO_g16_aave.120920.nc map_r05_TO_g16_aave.120920.nc
ln -sf ${USERINPUT}/cpl/cpl6/map_r05_to_gx1v6_e1000r300_090226.nc map_r05_to_gx1v6_e1000r300_090226.nc

cd ${MYINPUTDATA}/cpl/gridmaps
unlink fv1.9x2.5
ln -sf ${USERINPUT}/cpl/gridmaps/fv1.9x2.5 fv1.9x2.5
unlink gx1v6
ln -sf ${USERINPUT}/cpl/gridmaps/gx1v6 gx1v6



#2.4 pflotran input data (if not coupled run, commment out)
PETSC_DIR=/software/user_tools/current/cades-ccsi/petsc4pf/openmpi-1.10-gcc-5.3
PETSC_ARCH=arch-orcondo-openmpi-gcc53-debug
cd ${MYINPUTDATA}/pflotran/global
unlink pflotran_clm.in
ln -sf pflotran_clm_th-f19g16x15_v20160608.in pflotran_clm.in

cd ${MYINPUTDATA}/pflotran
unlink global
ln -sf ${USERINPUT}/pflotran/global global


#2.1 create case/run directory, if not yet
CLM_CASE_DIR=$CCSI_DIR/proj-shared/project_acme/cases
mkdir -p ${CLM_CASE_DIR}
CLM_RUN_DIR=${WORKDIR}/runs
mkdir -p ${CLM_RUN_DIR}

#3.2 Run python scripts to create/build/submit a run. Please refer to the python scripts 'runGRDCLM.py' with --help
cd ${CLM_COUPLED_MODEL}/clm4-pf-tools
python ./runGRDCLM.py \
                 --caseidprefix=clm0pf \
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
                 --hist_userdefined=clm_output_adspinup \
                 --np=320 --ppn=32 \
                 --walltime=02:00:00 \
                 --jobname=clm0ad --no_submit \
                 --compset=I1850CRUCLM45CN \
                 --res=f19_g16 \
                 --coldstart --ad_spinup --run_n=0100 --run_unit=date \
                 --surfdatafile=surfdata_1.9x2.5_simyr1850.nc \
                 --clm_pflotran \
                 --petsc_dir=${PETSC_DIR} \
                 --petsc_arch=${PETSC_ARCH} \
                 --pflotran_srcdir=${PFLOTRAN_COUPLED_MODEL} \
                 --pflotran_indir=global \
                 --clm_pf_colmode


## OPTION I: if coupled run with PFLOTRAN, uncomment the above 5 lines and add a '\' continuation sign properly
#                 --clm_pflotran \
#                 --petsc_dir=${PETSC_DIR} \
#                 --petsc_arch=${PETSC_ARCH} \
#                 --pflotran_srcdir=${PFLOTRAN_COUPLED_MODEL} \
#                 --pflotran_indir=global \
#                 --clm_pf_colmode \

## OPTION II: cold-start or having an 'initial case' for (ad)spinup, please ONLY having ONE of the following first 2 lines (modifying it properly)
#                 --finidat_case=clm1pfC40x40bgc_US-Brw_I1850CLM45CN_ad_spinup --finidat_year=1001 --run_n=2001 --run_unit=date \
#                 --coldstart --ad_spinup --run_n=1001 --run_unit=date \
#                 --compset=I1850CLM45CN \

# OPTION III: for transient run, OPTION III should be like the following
#                 --finidat_case=clm1pfC40x40bgc_US-Brw_I1850CLM45CN --finidat_year=2001 --run_n=2014 --run_unit=date \
#                 --compset=I20TRCLM45CN \
 

cd ${CASESETUP_DIR}

##----------------------------#
## this is for Global CLM4.5_r85 and clm45-pflotran model setup/runs on CADES-OR-CONDON
## Please direct your comments, questions or ask for instructions to
## 
## Fengming YUAN, CCSI/ESD-ORNL
## (865)57408486(O)
## yuanf@ornl.gov
## 
## 2017-02-25
##----------------------------#
