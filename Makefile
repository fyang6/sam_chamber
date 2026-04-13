# Makefile for various platforms
# Execute using Build csh-script only!
# Used together with Perl scripts in SRC/SCRIPT 
# (C) 2005 Marat Khairoutdinov
# $Id: Makefile 1595 2014-10-30 21:36:30Z dschanen@uwm.edu $
#------------------------------------------------------------------
# uncomment to disable timers:
#
#NOTIMERS=-DDISABLE_TIMERS
#-----------------------------------------------------------------

SAM = SAM_$(ADV_DIR)_$(SGS_DIR)_$(RAD_DIR)_$(MICRO_DIR)

# Determine platform 
PLATFORM := $(shell uname -s)

#------------------------------------------------------------------------
# This Makefile is for SAM running on Portage/Superior at MTU
# If having any question, contact Dr. Gowtham: g@mtu.edu
# MTU (portage, INTEL )
#
ifeq ($(PLATFORM),Linux)

NETCDF_DIR  = /home/fanyang/usr
MPICH_DIR   = /home/fanyang/lib/mpich
INC_NETCDF := $(NETCDF_DIR)/include
LIB_NETCDF := $(NETCDF_DIR)/lib
INC_MPI    := $(MPICH_DIR)/include
LIB_MPI    := $(MPICH_DIR)/lib

FF77        = mpif90 -c -fixed -extend_source
FF90        = mpif90 -c
CC          = mpicc -c -DLINUX

#FFLAGS      = -Os -pad
LD          = mpif90
FFLAGS      = -O3 
FFLAGS     += -r8
FFLAGS     += -I${INC_MPI} -I${INC_NETCDF} -DNETCDF
##CFFLAGS     = -I${INC_MPI} -I${INC_NETCDF} -DNETCDF
LDFLAGS     = -L${LIB_MPI} -lmpi -lmpifort -L${LIB_NETCDF} -lnetcdf -lnetcdff
##LDFLAGS    += -mkl

endif

#----------------------------------
#----------------------------------------------
# you don't need to edit below this line


#compute the search path
dirs := . $(shell cat Filepath)
VPATH    := $(foreach dir,$(dirs),$(wildcard $(dir))) 

.SUFFIXES:
.SUFFIXES: .f .f90 .F90 .F .c .o



all: $(SAM_DIR)/$(SAM)

SOURCES   := $(shell cat Srcfiles)

Depends: Srcfiles Filepath
	$(SAM_SRC)/SCRIPT/mkDepends Filepath Srcfiles > $@

Srcfiles: Filepath
	$(SAM_SRC)/SCRIPT/mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES))) 

$(SAM_DIR)/$(SAM): $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

.F90.o:
	${FF90}  ${FFLAGS} $(if $(filter $@, $(WARNABLE) ), $(WARNINGS) ) $<;
.f90.o:
	${FF90}  ${FFLAGS} $(if $(filter $@, $(WARNABLE) ), $(WARNINGS) ) $<;
.f.o:
	${FF77}  ${FFLAGS} $<
.F.o:
	${FF77}  ${FFLAGS} $<
.c.o:
	${CC}  ${CFLAGS} -I$(SAM_SRC)/TIMING $(NOTIMERS) $<

include Depends
