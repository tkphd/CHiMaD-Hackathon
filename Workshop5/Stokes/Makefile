# Makefile
# GNU makefile for example Stokes model with
# convex splitting and successive over-relaxation
# Questions/comments to trevor.keller@nist.gov (Trevor Keller)

# includes
incdir = $(MMSP_DIR)/include

# compilers/flags
compiler = g++
pcompiler = mpicxx
flags = -O3 -std=c++11 -I $(incdir)
pflags = $(flags) -include mpi.h
links = -lz -fopenmp

# the program
stokes: stokes.cpp
	$(compiler) $(flags) -fopenmp $< -o $@ $(links)

parallel: stokes.cpp
	$(pcompiler) $(pflags) $< -o $@ $(links)

clean:
	rm -f stokes parallel
