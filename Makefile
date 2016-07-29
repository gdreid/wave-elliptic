.IGNORE:

SHELL = /bin/sh

prefix = /usr/local
BBH_SYSTEM = LINUX_INTEL64
bindir = ../../bin

BBH_SYSTEM = LINUX_INTEL64

RNPL   = rnpl

AR     = xiar
RANLIB = ranlib

F77_TRANSFORM = touch
 
LIBS       = -lpamrw -lamrdw $(MPILIB) -lm -lmpi   -lbbhutil -lsvml -lifcore -lm 
LDFLAGS    = -L../../lib -O3    -L/usr/local/openmpi/lib -L/usr/X11R6/lib -L. -L$(prefix)/lib 


CC       = mpicc
CFLAGS     = -O3   
CPPFLAGS = -I. -I/usr/local/openmpi/include -I../../include -I$(prefix)/include \
             -I/usr/local/include

CC_COMP  = $(CC) $(CPPFLAGS) -c $(CFLAGS)
CC_LOAD  = $(CC) $(CFLAGS) $(LDFLAGS) 

F77      = mpif77
F77FLAGS =   

F77_COMP   = $(F77) -c $(F77FLAGS) 
F77_LOAD   = $(F77) $(F77FLAGS) $(F77_LDFLAGS) $(LDFLAGS) 

EXECUTABLES     = wave 
# Miscellaneous files to clean up
MISCDATAFILES   =

SRC = *.f *.inc

all: $(EXECUTABLES)

.f.o:
	$(F77_COMP) $*.f 

.c.o:
	$(CC_COMP) -c $*.c

all: $(EXECUTABLES)

install: all
	echo "Made all ... no installation"

full: install confidence_tests

confidence_tests:
	echo "Not implemented yet"

translate: 
	touch translate

wave.o: wave.c

WAVE_OBJS = wave.o init_f.o init_f_t.o res_f.o res_f_t.o u_f.o u_f_t.o res_psi_t0.o res_mg_psi_t0.o lop_psi_t0.o relax_psi_t0.o num.o
wave: $(WAVE_OBJS) 
	$(CC_LOAD) $(WAVE_OBJS) $(LIBS) -o wave

########################################################################
# Clean-up
########################################################################
clean:
	/bin/rm $(EXECUTABLES)
	/bin/rm *_.c > /dev/null 2>&1 
	/bin/rm *.o 
	/bin/rm *.a
	/bin/rm *.sdf
	/bin/rm *.segdat
	/bin/rm *~
	/bin/rm  config.cache config.log config.status
	/bin/rm $(MISCDATAFILES)
	/bin/rm run*/*.sdf
	/bin/rm run*/*.dat
	/bin/rm run*/mfile
