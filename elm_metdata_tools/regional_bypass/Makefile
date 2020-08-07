# ============================================================================
# Name        : Makefile
# Author      : 
# Version     :
# ============================================================================

#FC_ROOT := /usr/local/gcc-x/gcc-x-clang
#MPI_ROOT := /usr/local/gcc-x/mpich-3.3
#NETCDF_PATH := /usr/local/gcc-x/netcdf-4.x-hdf5

FC_ROOT := /software/dev_tools/swtree/cs400_centos7.5_pe2018/gcc/8.1.0/centos7.5_gnu4.8.5
MPI_ROOT := /software/dev_tools/swtree/cs400_centos7.5_pe2018/openmpi/3.1.5/centos7.5_gnu8.1.0
NETCDF_PATH := /software/user_tools/current/cades-ccsi/netcdf4/openmpi-3.1-gcc-8.1
PATH := $(FC_ROOT)/bin:$(PATH)

INCLUDE := -I$(NETCDF_PATH)/include -I$(MPI_ROOT)/include -I$(FC_ROOT)/include
LIBS := -L$(NETCDF_PATH)/lib -lnetcdff -lnetcdf -L$(MPI_ROOT)/lib -L$(FC_ROOT)/lib

PHONY: all clean_build

all: makezones_reanalysis.f90
	$(MPI_ROOT)/bin/mpif90 -O2 -g $(INCLUDE) $(LIBS) \
	    -o ./makezones \
	    ./makezones_reanalysis.f90

clean_build:
	rm -f *.mod *.o

clean: clean_build
	rm -f makezones
