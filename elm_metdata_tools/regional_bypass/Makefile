# ============================================================================
# Name        : Makefile
# Author      : 
# Version     :
# ============================================================================

FC_ROOT := /usr/local/gcc-clang-darwin/gcc-x-clang
MPI_ROOT := /usr/local/gcc-clang-darwin/openmpi-4.x-gcc
NETCDF_PATH := /usr/local/gcc-clang-darwin/netcdf-4.x-hdf5-openmpi
NETCDF_F_PATH := $(NETCDF_PATH)

#FC_ROOT := /software/dev_tools/swtree/cs400_centos7.5_pe2018/gcc/8.1.0/centos7.5_gnu4.8.5
#MPI_ROOT := /software/dev_tools/swtree/cs400_centos7.5_pe2018/openmpi/3.1.5/centos7.5_gnu8.1.0
#NETCDF_PATH := /software/user_tools/current/cades-ccsi/netcdf4/openmpi-3.1-gcc-8.1

#MPI_ROOT := /sw/andes/spack-envs/base/opt/linux-rhel8-x86_64/gcc-9.3.0/openmpi-4.0.4-skxqfeiocc5jtuw3y6dwtnxzzqjp5ffs
#FC_ROOT := /sw/andes/gcc/9.3.0
#NETCDF_PATH := /sw/andes/spack-envs/base/opt/linux-rhel8-x86_64/gcc-9.3.0/netcdf-c-4.7.4-dslyq2nupd5wb222cwn7p7qvka4ei5uk
#NETCDF_F_PATH := /sw/andes/spack-envs/base/opt/linux-rhel8-x86_64/gcc-9.3.0/netcdf-fortran-4.5.3-2ehk7dtlhazxhugbfmp77brtto3z2x5l

PATH := $(FC_ROOT)/bin:$(PATH)

INCLUDE := -I$(NETCDF_F_PATH)/include -I$(NETCDF_PATH)/include -I$(MPI_ROOT)/include -I$(FC_ROOT)/include
LIBS := -L$(NETCDF_F_PATH)/lib -lnetcdff -L$(NETCDF_PATH)/lib -lnetcdf -L$(MPI_ROOT)/lib -L$(FC_ROOT)/lib

PHONY: all clean_build

all: makezones_reanalysis.f90
	$(MPI_ROOT)/bin/mpif90 -O2 -g $(INCLUDE) $(LIBS) \
	    -o ./makezones \
	    ./makezones_reanalysis.f90

clean_build:
	rm -f *.mod *.o

clean: clean_build
	rm -f makezones
