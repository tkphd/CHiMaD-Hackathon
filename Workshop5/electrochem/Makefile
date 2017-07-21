# Makefile
# GNU makefile for example electrochem benchmark code
# Questions/comments to trevor.keller@nist.gov (Trevor Keller)

# includes
incdir = $(MMSP_PATH)/include

# compilers/flags
compiler = g++
flags = -O3 -Wall -std=c++11 -I $(incdir)
links = -lz -fopenmp
pcompiler = mpicxx
pflags = $(flags) -include mpi.h

# the program
electrochem: electrochem.cpp
	$(compiler) $(flags) $< -o $@ $(links)

parallel: electrochem.cpp
	$(pcompiler) $(pflags) $< -o $@ $(links)

clean:
	rm -f electrochem parallel
