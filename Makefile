# Makefile
# GNU makefile for example Cahn-Hilliard model code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = $(MMSP_PATH)/include
ANACONDA = /users/tnk10/.conda/envs/mmsp

# compilers/flags
compiler = g++
flags = -O3 -I $(incdir)

# dependencies
core = $(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp

# the program
mmsp2png: mmsp2png.cpp $(core) $(ANACONDA)/include/IL/devil_cpp_wrapper.hpp
	$(compiler) $(flags) -I $(ANACONDA)/include $< -o $@ -L $(ANACONDA)/lib -lz -lIL -lILU -lILUT

clean:
	rm -f mmsp2png
