# Makefile
# GNU makefile for example Cahn-Hilliard model code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = $(MMSP_PATH)/include
IL_PATH=/users/tnk10/.conda/envs/mmsp

# compilers/flags
compiler = g++
flags = -O3 -I $(incdir)

# dependencies
core = $(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp

# the program
mmsp2png: mmsp2png.cpp $(core) $(IL_PATH)/include/IL/devil_cpp_wrapper.hpp
	$(compiler) $(flags) -I$(IL_PATH)/include -include IL/il.h $< -o $@ -L$(IL_PATH)/lib -lz -lIL -lILU -lILUT

linescan: linescan.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

clean:
	rm -f mmsp2png linescan
