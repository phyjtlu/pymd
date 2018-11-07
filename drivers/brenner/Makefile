#
#.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90 .c

FOR = gfortran  
FFLAGS = -O3 -fno-automatic 
#FOR = ifort
#FFLAGS = -O3

#CC=icc
CC=gcc
LIBS = -lm

OBJS = main.o brenner_sub.o
MLOBJS = brenner_sub.o Brenner.o Brennertm.o
MLFILE = Brenner.c Brenner.tm

 



MLINKDIR = /applications/Mathematica/7.0/SystemFiles/Links/MathLink/DeveloperKit
SYS = Linux-x86-64
EXTRALIBS = -lm -lpthread -lrt # Set these with appropriate libs for your system.
MLLIB = ML64i3 # Set this to ML64i3 if using a 64-bit system

CADDSDIR = ${MLINKDIR}/${SYS}/CompilerAdditions
INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}


prep=${CADDSDIR}/mprep


# MathLink executable
Brenner : $(MLOBJS)
	$(CC) -I${INCDIR} $(MLOBJS) -L${LIBDIR} -l${MLLIB} ${EXTRALIBS} $(LIBS) -o $@

# Normal Linux executable
BrennerMD : $(OBJS) 
	$(FOR) $(OBJS) $(LIBS) -o BrennerMD
# Matlab executable
brennermex : brennermex.F
	mex brennermex.F
# Python lib
BrennerPY : brenner_sub.f
	f2py --fcompiler=$(FOR) -c brenner_sub.f -m brenner

#
Brennertm.c : Brenner.tm
	        ${prep} Brenner.tm -o Brennertm.c

.c.o:
	        ${CC} -c -I${INCDIR} $<
.f.o:
		$(FOR) $(FFLAGS) -c $<

clean :
	rm -f *.o Brenner BrennerMD brennermex.mexa64 brenner.so
