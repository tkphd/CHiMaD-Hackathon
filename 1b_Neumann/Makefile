# Makefile
# GNU makefile for example Cahn-Hilliard model code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = $(MMSP_PATH)/include

# compilers/flags
compiler = g++
flags = -O3 -I $(incdir)
pcompiler = mpicxx
ccompiler = mpicxx.mpich2
pflags = -O3 -I $(incdir) -include mpi.h

# dependencies
core = $(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp

# the program
cahn-hilliard: cahn-hilliard.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

parallel: cahn-hilliard.cpp $(core)
	$(pcompiler) $(pflags) $< -o $@ -lz

ctcms: cahn-hilliard.cpp $(core)
	$(ccompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f cahn-hilliard parallel
