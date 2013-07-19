########################################################################
# Compiler and external dependences
########################################################################
CC        = icc
F77       = ifort
CXX       = icpc
F90       = ifort
#HYPRE_DIR = /home3/idobbsdi/hypre/hypre-2.6.0b-babel/src/hypre
#HYPRE_DIR = /home3/idobbsdi/hypre/hypre-2.4.0b/src/hypre
#HYPRE_DIR = /home3/idobbsdi/hypre/hypre-2.2.0b/src/hypre
HYPRE_DIR = /home3/idobbsdi/hypre/h2.2_withnewfiles/hypre-2.2.0b/src/hypre
#HYPRE_DIR = /nasa/hypre/2.2.0b/intel
#HYPRE_DIR = /nasa/hypre/2.2.0b/intel-10.1.013_64/mpt.1.23.pre

########################################################################
# Compiling and linking options
########################################################################
COPTS     = -g -pedantic -Wall
CINCLUDES = -I$(HYPRE_DIR)/include
CDEFS     = -DHAVE_CONFIG_H -DHYPRE_TIMING
CFLAGS    = $(COPTS) $(CINCLUDES) $(CDEFS)
#FOPTS     = -g -O2
FOPTS     = -O2
#FOPTS     = -O3
#FOPTS     = -O3 -ipo
FINCLUDES = $(CINCLUDES)
FDEFS     = $(CDEFS)
FFLAGS    = $(FOPTS) $(FINCLUDES) $(FDEFS)
CXXOPTS   = $(COPTS)
CXXINCLUDES = $(CINCLUDES) -I..
CXXDEFS   = $(CDEFS)
#IFLAGS_BXX = -I$(HYPRE_DIR)/babel-runtime/sidl
#CXXFLAGS  = $(CXXOPTS) $(CXXINCLUDES) $(CXXDEFS) $(IFLAGS_BXX)
#IF90FLAGS = -I$(HYPRE_DIR)/babel/bHYPREClient-F90
#F90FLAGS = $(FFLAGS) $(IF90FLAGS)
LIBS      = -L$(HYPRE_DIR)/lib -lHYPRE -lg2c -lm
LFLAGS    = $(LINKOPTS) $(LIBS) -lstdc++


#--this works
LFLAGS_B = -L${HYPRE_DIR}/lib -lbHYPREClient-C -lbHYPREClient-CX -lbHYPREClient-F -lbHYPRE -lsidl -ldl -lHYPRE -lm -lstdc++ -lmpi
#-lxml2
#The below flags are helpful for looking for floating-point assist faults
#TRACEBACKFLAG = -O0 -fpe0 -traceback -g
#Debugging flags (from nasa web seminar)
#DEBUGFLAG = -g -traceback -check -fpe0
#DEBUGFLAG = -g


########################################################################
# Rules for compiling the source files
########################################################################
.SUFFIXES: .c .f

.c.o:
	$(CC) $(CFLAGS) -c $<
.f.o:
	$(F77) $(FFLAGS) $(TRACEBACKFLAG) $(DEBUGFLAG) -c $<

########################################################################
# TRIAL PROGRAM
########################################################################

FILENAMES= ppmod.o ppinit.o pp.o ppbou.o ppphy.o ppsub.o ppfor.o \
	ppadv.o ppprt.o pprad.o ppmat.o pprphy.o \
	ppvisc.o ppshift.o ppenergy.o ppsor.o pptest.o \
	ppcool.o ppscalar.o ppdif.o ppirad.o ppsaum.o ppplanets.o pphypre.o \
	ppchemistry.o ppheld.o ppeliza.o ppwavelth.o
# laplacetest.o 
default::start

start: $(FILENAMES)
	$(F77) -o ttx $(FILENAMES) $(TRACEBACKFLAG) $(DEBUGFLAG) $(LFLAGS_B) -lmpi

########################################################################
# Clean up
########################################################################
clean:
	rm -f *.mod
	rm -f *.o
	rm -f *~
