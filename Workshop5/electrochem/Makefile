# Makefile
# GNU makefile for example electrochem benchmark code
# Questions/comments to trevor.keller@nist.gov (Trevor Keller)

# includes
incdir = $(MMSP_PATH)/include

# compilers/flags
compiler = g++
flags = -O3 -Wall -pedantic -std=c++11 -I $(incdir)
links = -lz -fopenmp
pcompiler = mpicxx
pflags = $(flags) -include mpi.h

# the program
electrochem: electrochem.cpp electrochem.hpp energy.hpp
	$(compiler) $(flags) $< -o $@ $(links)

ichem: electrochem.cpp electrochem.hpp energy.hpp
	icc $(flags) -w3 $< -o $@ $(links)

parallel: electrochem.cpp electrochem.hpp energy.hpp
	$(pcompiler) $(pflags) $< -o $@ $(links)

# utilities
linescan: linescan.cpp
	$(compiler) -O2 -Wall -pedantic -std=c++11 -I $(incdir) $< -o $@ -lz

clean:
	rm -f electrochem ichem linescan parallel
