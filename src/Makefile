# Compiler:
F90C  = mpif90
CC  = mpicc
F90=ifort

# Link libraries and options
FFLAGS  = -O3 -cm -w -Wall -vec_report0 -sox -openmp -openmp_report0 -I$(HEALPIX)/include 
FFLAGS_PAR = -03 -heap-arrays 16 -I$(HEALPIX)/include  -I$(CFITSIO)/include 
CFLAGS	=-O2 -fomit-frame-pointer -fforce-addr -fforce-mem
LDFLAGS	= -L$(HEALPIX)/lib -L$(CFITSIO)/lib -lhealpix -lgif -lcfitsio 


#LDFLAGS = -lmpi -lfmpi 
#INTLIB = -L/mn/stornext/u2/yabebalf/lib -lintlib
#FILLIB =  -L/mn/stornext/u2/yabebalf/lib -lfilonsimp


# Set up appropriate suffix list
.SUFFIXES:
.SUFFIXES: .o .f .for .f90 .c

# Compile-time options
# F90 options

BIN=../bin


# Obj list:
MOD	=  healpix_types.o parsers.o 
OBJ	=  box2ball.o

MODT	=  parsers.o common_module.o
OBJT	=  test_filon_quad.o

MODT1	=  parsers.o common_module.o
OBJT1	=  test_fil_simp_points.o 

MOD1 = parsers.o common_module.o
OBJ1 = ball2pix.o

#+
# Make some default rules
#-

# Double suffix rules:
.f90.o .f.o .for.o :
	$(F90C) $(FFLAGS) -c $<
.c.o	:
	$(CC) $(CFLAGS) -c $<

# Single suffix rules:
.f90 .f .for :
	$(F90C) $(FFLAGS) -o $@ $<

.o :
	$(F90C) $(LDFLAGS) -o $@ $<



# Build all programs
all: box2ball ball2pix
test: test_filon_quad test_filsimp

test_filsimp: ${MODT1} ${OBJT1}
	$(F90C) $(FFLAGS) -o ${BIN}/$@ ${OBJT1} ${MODT1} ${LDFLAGS}

test_filon_quad: ${MODT} ${OBJT}
	$(F90C) $(FFLAGS) -o ${BIN}/$@ $@.f90 ${MODT} ${LDFLAGS} 

box2ball: ${MOD} ${OBJ}
	$(F90C) $(FFLAGS) -o ${BIN}/$@ ${MOD} ${OBJ} 

ball2pix: ${MOD1}
	$(F90C) $(FFLAGS) -o ${BIN}/$@ $@.f90 ${MOD1} ${LDFLAGS} 


# Remove all files created from building 
clean	:
	-$(RM) hess
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.lst

tidy:	clean

#+
# Dependencies
#-


parsers.o : parsers.f90
healpix_types.o : healpix_types.f90
box2ball.o : box2ball.f90
ball2pix.o : ball2pix.f90