# Compiler:
F90C  = mpif90
CC  = mpicc
F90=ifort

# Link libraries and options
FFLAGS  = -O3 -cm -w -Wall -vec_report0 -sox -openmp -openmp_report0 -I$(HEALPIX)/include 
FFLAGS_PAR = -03 -heap-arrays 16 -I$(HEALPIX)/include  -I$(CFITSIO)/include 
CFLAGS	=-O2 -fomit-frame-pointer -fforce-addr -fforce-mem
LDFLAGS	= -L$(HEALPIX)/lib -L$(CFITSIO)/lib -lhealpix -lgif -lcfitsio 

ifeq ($(TARGET),nersc) 
F90C =  ftn #mpif90	
F90 = ftn #mpif90

# Compile-time options
HEALPIX=$(HEALPIX3)
# F90 options  ## -O3 -xAVX -mavx 
#check http://www.nas.nasa.gov/hecc/support/kb/recommended-intel-compiler-debugging-options_92.html
DEBUG = -g -debug inline-debug-info -ftrapuv 
RUNINFO=-check bounds -traceback -traceback -gen-interfaces  #-warn none 
FFLAGS = $(DEBUG) $(RUNINFO) -heap-arrays -I$(HEALPIX)/include  -I$(CFITSIO)/include  -qopenmp
# CC options:
CFLAGS	= #-O2 -fomit-frame-pointer -fforce-addr -fforce-mem
# Link libraries and options
LDFLAGS	= -L$(HEALPIX)/lib -L$(CFITSIO)/lib -lhealpix -lgif -lcfitsio 
LAPACK  = -I$(MKLROOT)/include -mkl
LAPACKLIB = $(LAPACK) 
endif


#LDFLAGS = -lmpi -lfmpi 
#INTLIB = -L/mn/stornext/u2/yabebalf/lib -lintlib
#FILLIB =  -L/mn/stornext/u2/yabebalf/lib -lfilonsimp


# Set up appropriate suffix list
.SUFFIXES:
.SUFFIXES: .o .f .for .f90 .c

# Compile-time options
# F90 options

BIN=../bin

MOD_C1=parsers.o mapio.o
MOD_C2=fftw_wrapper.o almnutil.o needlet_kernel.o almn_tools.o
MOD_C12=$(MOD_C1) $(MOD_C2)

# Obj list:
MOD_BOX2BALL=$(MOD_C1)
OBJ_BOX2BALL=box2ball.o

MOD_BALL2PIX=$(MOD_C1)
OBJ_BALL2PIX=ball2pix.o

#ball2almn and almn2ball cases
MOD_BALL2ALMN=$(MOD_C12)
OBJ_BALL2ALMN=ball2almn.o

MOD_ALMN2BALL=$(MOD_C12)
OBJ_ALMN2BALL=almn2ball.o

#almn2beta and beta2almn cases
MOD_ALMN2BETA=$(MOD_C12)
OBJ_ALMN2BETA=almn2beta.o

#ball2beta and beta2ball cases
MOD_BALL2BETA=$(MOD_C12)
OBJ_BALL2BETA=ball2beta.o

MOD_BETA2BALL=$(MOD_C12)
OBJ_BETA2BALL=beta2ball.o

# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(F90C) $(FFLAGS) -c $<

%.o: %.F90
	$(F90C) $(FFLAGS) -c $<

# Make some default rules
# Double suffix rules:
.f.o .for.o :
	$(F90C) $(FFLAGS) -c $<
.c.o	:
	$(CC) $(CFLAGS) -c $<

# Single suffix rules:
.f90 .F90 .f .for :
	$(F90C) $(FFLAGS) -o $@ $<

.o :
	$(F90C) $(LDFLAGS) -o $@ $<



# Build all programs
default: ball2almn almn2ball almn2beta beta2ball
all: box2ball ball2pix almr2cl $(default)

box2ball: ${MOD_BOX2BALL} ${OBJ_BOX2BALL}
	$(F90C) $(FFLAGS) -o ${BIN}/$@ $@.f90 ${MOD_BOX2BALL} ${LDFLAGS} 

ball2pix: ${MOD_BALL2PIX} ${OBJ_BALL2PIX}
	$(F90C) $(FFLAGS) -o ${BIN}/$@ $@.f90 ${MOD_BALL2PIX} ${LDFLAGS} 

ball2almn: ${MOD_BALL2ALMN} ${OBJ_BALL2ALMN}
	$(F90C) $(FFLAGS) -o ${BIN}/$@ $@.f90 ${MOD_BALL2ALMN} ${LDFLAGS} 

almn2ball: ${MOD_ALMN2BALL} ${OBJ_ALMN2BALL}
	$(F90C) $(FFLAGS) -o ${BIN}/$@ $@.f90 ${MOD_ALMN2BALL} ${LDFLAGS} 

almn2beta: ${MOD_ALMN2BETA} ${OBJ_ALMN2BETA}
	$(F90C) $(FFLAGS) -o ${BIN}/$@ $@.f90 ${MOD_ALMN2BETA} ${LDFLAGS} 

ball2beta: ${MOD_BALL2BETA} ${OBJ_BALL2BETA}
	$(F90C) $(FFLAGS) -o ${BIN}/$@ $@.f90 ${MOD_BALL2BETA} ${LDFLAGS} 

beta2ball: ${MOD_BETA2BALL} ${OBJ_BETA2BALL}
	$(F90C) $(FFLAGS) -o ${BIN}/$@ $@.f90 ${MOD_BETA2BALL} ${LDFLAGS} 

almr2cln:  ${MOD_C12}
	$(F90C) $(FFLAGS) -o ${BIN}/$@ $@.f90 ${MOD_C12} ${LDFLAGS} 

idlball:
	echo "compiling idlball directory"
	$(MAKE) -C $@


#++++++++++++++++++++++++++++++++++++++++++++++++++
# Dependencies
#--------------------------------------------------


almn_tools.o : almn_tools.F90
parsers.o : parsers.f90
healpix_types.o : healpix_types.f90
box2ball.o : box2ball.f90
ball2pix.o : ball2pix.f90
ball2almn.o: ball2almn.f90
almn2ball.o: almn2ball.f90
ball2beta.o: ball2beta.f90
beta2ball.o: beta2ball.f90

# Utility targets
.PHONY: clean veryclean

veryclean: clean
	rm -f *~ ~# ${BIN}/*

# Remove all files created from building 
clean	:
	-$(RM) hess
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.lst

tidy:	clean
