#by default parallel, use 'make SERIAL=1' for serial

ifneq ($(SERIAL),1)
CC=mpicc 
CFLAGS=-O3 -fopenmp -DMPI -I$(TACC_GSL_INC) -I$(TACC_GSL_INC)/gsl -I$(TACC_HDF5_INC) 

else
CC=h5cc
CFLAGS=-O2 -DNOSILO -fopenmp -I$(TACC_GSL_INC) -I$(TACC_GSL_INC)/gsl -I$(TACC_HDF5_INC) #-DMPI

//CC=clang
//CFLAGS = -O2 -Wno-unused-result -I/usr/lib/gcc/x86_64-linux-gnu/5.4.0/include -I/usr/include/hdf5/serial -Wunused-function -fopenmp=libiomp5 -g 
-fsanitize=address -fno-omit-frame-pointer

//CC=/usr/bin/h5cc
//CFLAGS = -O2 -Wno-unused-result -I/usr/lib/gcc/x86_64-linux-gnu/5.4.0/include -I/usr/include/hdf5/serial -Wunused-function -fopenmp 

endif

LIBS=-lm -lgsl -lgslcblas -lrt -L$(TACC_GSL_LIB) -lhdf5 -L$(TACC_HDF5_LIB)
RM=/bin/rm
OBJS = mpi.o u2prad.o magn.o postproc.o fileop.o misc.o physics.o finite.o problem.o metric.o relele.o rad.o opacities.o u2p.o frames.o p2u.o nonthermal.o u2p_ff.o 

all: ko phisli thsli phiavg regrid ana avg

ko: ko.o $(OBJS) Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o ko ko.o $(OBJS) $(LIBS)

ana: ana.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o ana ana.o $(OBJS) $(LIBS)

avg: avg.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o avg avg.o $(OBJS) $(LIBS)

outavg: outavg.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o outavg outavg.o $(OBJS) $(LIBS)

phisli: phisli.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o phisli phisli.o $(OBJS) $(LIBS)

thsli: thsli.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o thsli thsli.o $(OBJS) $(LIBS)

phiavg: phiavg.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o phiavg phiavg.o $(OBJS) $(LIBS)

regrid: regrid.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o regrid regrid.o $(OBJS) $(LIBS)

#dumps2hdf5: dumps2hdf5.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
#	$(CC) $(CFLAGS) -o dumps2hdf5 dumps2hdf5.o $(OBJS) $(LIBS)

clean:
	$(RM) -f ko ana avg phiavg phisli thsli outavg regrid *~ *.o *.oo
