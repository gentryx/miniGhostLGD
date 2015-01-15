#@header
# ************************************************************************
#
#      miniGhost: stencil computations with boundary exchange.
#              Copyright (2012) sandia corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Richard F. Barrett (rfbarre@sandia.gov) or
#                    Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
#@header

# Simple hand-tuned makefile, configured for the gnu compiler and MPICH.
# Modify as necessary for your environment.

#-----------------------------------------------------------------------

# Option: -D_MG_MPI
PROTOCOL = -D_MG_MPI
#PROTOCOL = -D_MG_SERIAL
#PROTOCOL += -D_MG_LGD
# Compilers
FC=mpif90
#FC=gfortran
#FC=ftn
CC=mpicc
#CC=cc
CXX=g++

CFLAGS = $(PROTOCOL)
# C main calling Fortran subroutine:
CFLAGS += -Df2c_

FFLAGS = $(PROTOCOL)

# Variable precision: -D_INT8 and/or -D_REAL8.
FFLAGS += -D_MG_INT4 -D_MG_REAL8
# Optimization
# Next line PGI compiler
#OPT_F = -fast -fastsse
#OPT_F = -g
#FFLAGS += $(OPT_F)

# Free-form Fortran source code:
# Next line PGI compiler:
#FFLAGS += -Mfree
# Next line Intel compiler:
#FFLAGS += -free -w
# Next line Cray CCE compiler:
#FFLAGS += -f free -m 4
# Next line gfortran
FFLAGS += -ffree-form

#disables line length restrictions which prevent us from compiling 
FFLAGS += -ffree-line-length-0
#FFLAGS += -cpp

CPPFLAGS = -std=c++11  -Wall -Wno-sign-promo -Wnon-virtual-dtor -march=native -fopenmp -O3 -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -I/usr/local/include/libflatarray/ -I/usr/local/include/libgeodecomp/ -I/usr/include/qt4/ 

# Array bounds checking: (expensive!)
#FFLAGS += -Mbounds
# Compile to include checkpointing capability.
#FFLAGS += -D_MG_CHECKPT

LD = $(CXX)
LDFLAGS = $(CFLAGS)
LIBS =

LIBS += -rdynamic 
#-----QTCore ?
LIBS += -lgeodecomp
LIBS += -lboost_date_time -lboost_filesystem -lboost_system
LIBS += -lboost_mpi
LIBS += -lpthread
LIBS += -lmpi_cxx -lmpi_f90 -lmpi_f77 -lmpi
LIBS += -ldl
LIBS += -lhwloc
LIBS += -lgfortran 
#LIBS += -lstdc++
LIBS += -lquadmath
LIBS += -Wl,-rpath,/usr/local/lib
#LIBS += -Wl,-rpath,/usr/lib/i386-linux-gnu -Wl,-rpath,/usr/lib -L/usr/local/lib

# g++ -o miniGhost.x main.o wrapper.o MG_CONSTANTS.o MG_OPTIONS.o MG_UTILS.o MG_ENV.o MG_PROFILING.o MG_GET_FACE.o MG_BUFINIT.o MG_PACK.o MG_IRECV.o MG_UNPACK_SVAF.o MG_SEND_BSPMA.o MG_SEND_SVAF.o MG_UNPACK_BSPMA.o MG_BSPMA.o MG_SVAF.o MG_FLUX_ACCUMULATE.o MG_STENCIL_COMPS.o MG_STENCIL.o MG_BSPMA_DIAGS.o MG_SVAF_DIAGS.o MG_ALLREDUCE.o MG_SUM_GRID.o MG_CHECKPOINT.o DRIVER_BSPMA.o DRIVER_SVAF.o DRIVER.o -rdynamic /usr/lib/i386-linux-gnu/libQtCore.so /usr/local/lib/libgeodecomp.so -lboost_date_time -lboost_filesystem -lboost_system -lboost_mpi -lpthread -lmpi_cxx -lmpi_f90 -lmpi_f77 -lmpi -ldl -lhwloc -lgfortran -lquadmath -Wl,-rpath,/usr/local/lib

include make_targets

# End makefile

