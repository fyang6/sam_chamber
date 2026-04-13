# Makefile to compile conversion utilities. You need to set include and library paths for NETCDF
#

PLATFORM := $(shell uname -s)

# uncomment line for your platform or add your own:

#---------------------------------------------------------
# Cray 

ifeq ($(PLATFORM),Linux)
NETCDF_DIR = /home/fanyang/usr
INC_NETCDF := $(NETCDF_DIR)/include
LIB_NETCDF := $(NETCDF_DIR)/lib

FF = ifort -O3 -fixed -extend_source -I${INC_NETCDF}
LD = ifort

LDFLAGS = -L${LIB_NETCDF} -lnetcdf -lnetcdff
endif
#---------------------------------------------------------
#
# Linux Cluster (seawulf@sunysb)
#

#INC_NETCDF = /nfs/user04/marat/local/include
#LIB_NETCDF = /nfs/user04/marat/local/lib


#FF = ifort -O3 -fixed -extend_source -I${INC_NETCDF}
#LDFLAGS = -L${LIB_NETCDF} -lnetcdf
#---------------------------------------------------------

VPATH = ./SRC

all: bin2D2nc bin3D2nc bin3D2nc_micro stat2nc 

.f:   
	$(FF) -o $@ -I./SRC $< ./SRC/hbuf_lib.f ./SRC/cape.f ./SRC/cin.f $(LDFLAGS) 

clean: 
	rm bin* com* stat* 2* isccp* modis* misr* *.o glue*
