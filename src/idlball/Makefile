# -*- Makefile -*-
# ------------------------------ variables for F90 --------------------------

# Compiler:
F90C  = mpif90
CC  = mpicc
F90=ifort


# Link libraries and options
FFLAGS  = -O3 -cm -w -Wall -vec_report0 -sox -openmp -openmp_report0 -I$(HEALPIX)/include 
FFLAGS_PAR = -03 -heap-arrays 16 -I$(HEALPIX)/include  -I$(CFITSIO)/include 
CFLAGS	=-O2 -fomit-frame-pointer -fforce-addr -fforce-mem
LDFLAGS	= -L$(HEALPIX)/lib -L$(CFITSIO)/lib -lhealpix -lgif -lcfitsio 

#--------------------- case NERSC -------
ifeq ($(TARGET),nersc) 
F90C =  ftn #mpif90	
F90 = ftn #mpif90

# Compile-time options
HEALPIX=$(HEALPIX3)
# F90 options  ## -O3 -xAVX -mavx -heap-arrays
FFLAGS =  -openmp -I$(HEALPIX)/include  -I$(CFITSIO)/include 
# CC options:
CFLAGS	= #-O2 -fomit-frame-pointer -fforce-addr -fforce-mem
# Link libraries and options
LDFLAGS	= -L$(HEALPIX)/lib -L$(CFITSIO)/lib -lhealpix -lgif -lcfitsio 
LAPACK  = -I$(MKLROOT)/include -mkl
LAPACKLIB = $(LAPACK) 
endif


FC=$(F90C)


# Set up appropriate suffix list
.SUFFIXES:
.SUFFIXES: .o .f .for .f90 .c

RM	= rm -f


# ------------------------------ variables --------------------------

#VPATH	= $(HEALPIX)/src/idl/almball
VPATH=../../pro


OBJ=../common_module.f90

exectest	= $(VPATH)/test_filon_quad
srctest	= test_filon_quad.f90
objtest 	= wav_ball_mod.o test_filon_quad.o fft.o 
inctest     = 

#ball2almn and almn2ball cases
exec1	= $(VPATH)/idl_almn2ball
src1	= idl_almn2ball.f90
obj1 	= wav_ball_mod.o idl_almn2ball.o fft.o
inc1     = 

exec2	= $(VPATH)/idl_ball2almn
src2	= idl_ball2almn.f90
obj2 	= wav_ball_mod.o idl_ball2almn.o fft.o
inc2     = 


#almn2beta and beta2almn cases
exec3	= $(VPATH)/idl_almn2beta
src3	= idl_almn2beta.f90
obj3 	= wav_ball_mod.o idl_almn2beta.o fft.o
inc3    = 

exec3b	= $(VPATH)/idl_almn2beta_par
src3b	= idl_almn2beta_par.f90
obj3b	= wav_ball_mod.o idl_almn2beta_par.o fft.o
inc3b	=

#idl_ball2beta_par and idl_beta2ball_par cases
exec4c	= $(VPATH)/idl_ball2beta_par
src4c	= idl_ball2beta_par.f90
obj4c 	= wav_ball_mod.o idl_ball2beta_par.o fft.o
inc4c    = 

exec4b	= $(VPATH)/idl_beta2ball_par
src4b	= idl_beta2ball_par.f90
obj4b 	= wav_ball_mod.o idl_beta2ball_par.o fft.o
inc4b    = 

#speed test healpix and seedlet
execS1	= $(VPATH)/speed_test_bl2bt_bt2bl
srcS1	= speed_test_bl2bt_bt2bl.f90
objS1 	= wav_ball_mod.o speed_test_bl2bt_bt2bl.o fft.o
incS1     = 

execS2	= $(VPATH)/speed_test_b2a_a2b
srcS2	= speed_test_b2a_a2b.f90
objS2 	= wav_ball_mod.o speed_test_b2a_a2b.o fft.o
incS2     = 

execS3	= $(VPATH)/speed_almn2ball_ball2almn
srcS3	= speed_almn2ball_ball2almn.f90
objS3 	= wav_ball_mod.o speed_almn2ball_ball2almn.o fft.o
incS3     = 

execS4	= $(VPATH)/speed_healpix
srcS4	= speed_healpix.f90
objS4 	= wav_ball_mod.o speed_healpix.o fft.o
incS4     = 


#sin cos basis for almn2ball and ball2almn cases
exec1c	= $(VPATH)/idl_almn2ball_sincos
src1c	= idl_almn2ball_sincos.f90
obj1c	= wav_ball_mod.o idl_almn2ball_sincos.o fft.o
inc1c   = 

exec2c	= $(VPATH)/idl_ball2almn_sincos
src2c	= idl_ball2almn_sincos.f90
obj2c 	= wav_ball_mod.o idl_ball2almn_sincos.o fft.o
inc2c     = 
#sin cos basis for ball2beta and beta2ball cases
exec4	= $(VPATH)/idl_ball2beta_sincos
src4	= idl_ball2beta_sincos.f90
obj4 	= wav_ball_mod.o idl_ball2beta_sincos.o fft.o
inc4    = 

exec4bb	= $(VPATH)/idl_beta2ball_sincos
src4bb	= idl_beta2ball_sincos.f90
obj4bb 	= wav_ball_mod.o idl_beta2ball_sincos.o fft.o
inc4bb    = 


#1-dimension cases
exec1b	= $(VPATH)/idl_an2ball
src1b	= idl_an2ball.f90
obj1b 	= wav_ball_mod.o idl_an2ball.o fft.o
inc1b     = 

exec2b	= $(VPATH)/idl_ball2an
src2b	= idl_ball2an.f90
obj2b 	= wav_ball_mod.o idl_ball2an.o fft.o
inc2b     = 


#other cases 
exec44c	= $(VPATH)/idl_psi_jqk
src44c	= idl_psi_jqk.f90
obj44c 	= wav_ball_mod.o idl_psi_jqk.o fft.o
inc44c    = 

exec4d	= $(VPATH)/idl_get_3dneedlet
src4d	= idl_get_3dneedlet.f90
obj4d 	= wav_ball_mod.o idl_get_3dneedlet.o fft.o
inc4d    = 

exec4cc	= $(VPATH)/idl_ball2ball
src4cc	= idl_ball2ball.f90
obj4cc 	= wav_ball_mod.o idl_ball2ball.o fft.o
inc4cc    = 

exec5	= $(VPATH)/idl_get_seedlet
src5	= idl_get_seedlet.f90
obj5 	= wav_ball_mod.o idl_get_seedlet.o fft.o
inc5    = 

exec6	= $(VPATH)/idl_beta2almn
src6	= idl_beta2almn.f90
obj6 	= idl_beta2almn.o
inc6    = 

exec6b	= $(VPATH)/idl_beta2almn_par
src6b	= idl_beta2almn_par.f90
obj6b 	= idl_beta2almn_par.o
inc6b    = 


# Double suffix rules:
.f90.o .f.o .for.o :
	$(F90C) $(FFLAGS) -c $(OBJ) $<
.c.o	:
	$(CC) $(CFLAGS) -c $(OBJ) $<

# Single suffix rules:
.f90 .f .for :
	$(F90C) $(FFLAGS) -o $@ $<

.o :
	$(F90C) $(LDFLAGS) -o $@ $<

#these are for sin and cosine routines

idl_par:  $(exec1) $(exec2) idl_beta2ball_par #idl_ball2ball_par
speed: speed_test_a2b_b2a speed_test_bl2bt_bt2bl
default: $(exec1c) $(exec2c) $(exec1b) $(exec2b)  ball2beta_sincos beta2ball_sincos get_3dneedlet psi_jqk


$(exectest): $(OBJ) $(objtest) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exectest) $(OBJ) $(objtest) $(LDFLAGS)


$(execS1): $(OBJ) $(objS1) #$(LIBDIR)/libhealpix.a
	ifort $(FFLAGS) -o $(execS1) $(OBJ) $(objS1) $(LDFLAGS)

$(execS2): $(OBJ) $(objS2) #$(LIBDIR)/libhealpix.a
	ifort $(FFLAGS) -o $(execS2) $(OBJ) $(objS2) $(LDFLAGS)

$(execS3): $(OBJ) $(objS3) #$(LIBDIR)/libhealpix.a
	ifort $(FFLAGS) -o $(execS3) $(OBJ) $(objS3) $(LDFLAGS)

$(execS4): $(OBJ) $(objS4) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(execS4) $(OBJ) $(objS4) $(LDFLAGS)



$(exec1): $(OBJ) $(obj1) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec1) $(OBJ) $(obj1) $(LDFLAGS)

$(exec2): $(OBJ) $(obj2) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec2) $(OBJ) $(obj2) $(LDFLAGS)

$(exec1c): $(OBJ) $(obj1c) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec1c) $(OBJ) $(obj1c) $(LDFLAGS)

$(exec2c): $(OBJ) $(obj2c) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec2c) $(OBJ) $(obj2c) $(LDFLAGS)

$(exec1b): $(OBJ) $(obj1b) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec1b) $(OBJ) $(obj1b) $(LDFLAGS)

$(exec2b): $(OBJ) $(obj2b) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec2b) $(OBJ) $(obj2b) $(LDFLAGS)

$(exec3): $(OBJ) $(obj3) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec3) $(OBJ) $(obj3) $(LDFLAGS)

$(exec3b): $(OBJ) $(obj3b) 
	$(FC) $(FFLAGS) -o $(exec3b) $(OBJ) $(obj3b) $(LDFLAGS)

$(exec4): $(OBJ) $(obj4) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec4) $(OBJ) $(obj4) $(LDFLAGS)

$(exec4b): $(OBJ) $(obj4b) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec4b) $(OBJ) $(obj4b) $(LDFLAGS)

$(exec4bb): $(OBJ) $(obj4bb) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec4bb) $(OBJ) $(obj4bb) $(LDFLAGS)

$(exec4c): $(OBJ) $(obj4c) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec4c) $(OBJ) $(obj4c) $(LDFLAGS)

$(exec44c): $(OBJ) $(obj44c) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec44c) $(OBJ) $(obj44c) $(LDFLAGS)

$(exec4d): $(OBJ) $(obj4d) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec4d) $(OBJ) $(obj4d) $(LDFLAGS)

$(exec4cc): $(OBJ) $(obj4cc) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec4cc) $(OBJ) $(obj4cc) $(LDFLAGS)

$(exec5): $(OBJ) $(obj5) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec5) $(OBJ) $(obj5) $(LDFLAGS)

$(exec6): $(OBJ) $(obj6) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec6) $(OBJ) $(obj6) $(LDFLAGS)

$(exec6b): $(OBJ) $(obj6b) #$(LIBDIR)/libhealpix.a
	$(FC) $(FFLAGS) -o $(exec6b) $(OBJ) $(obj6b) $(LDFLAGS)

clean:
	$(RM) *.o *.$(F90_MOD) lib*.a *.pc *.pcl *.il 
