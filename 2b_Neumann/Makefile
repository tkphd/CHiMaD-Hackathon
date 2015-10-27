# Makefile
# GNU makefile for example Ostwald ripening code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = $(MMSP_PATH)/include

# compilers/flags
compiler = g++
flags = -O3 -Wall -I $(incdir)
gflags = -g -pg -Wall -O2 -I $(incdir)
pcompiler = mpic++
ccompiler = mpic++.mpich2
fcompiler = /usr/mpi/gcc/openmpi-1.8.4rc1/bin/mpic++
pflags = -O3 -Wall -I $(incdir) -include mpi.h
fflags = -O3 -I $(incdir) -include /usr/mpi/gcc/openmpi-1.8.4rc1/include/mpi.h

# dependencies
core = $(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp

# the program
ostwald: ostwald.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

debug: ostwald.cpp $(core)
	$(compiler) -DDEBUG $(gflags) $< -o $@ -lz

parallel: ostwald.cpp $(core)
	$(pcompiler) $(pflags) $< -o $@ -lz

ctcms: ostwald.cpp $(core)
	$(ccompiler) $(pflags) $< -o $@ -lz

fnbnd: ostwald.cpp $(core)
	$(fcompiler) $(fflags) $< -o $@ -lz

clean:
	rm -f ostwald parallel ctcms debug fnbnd
