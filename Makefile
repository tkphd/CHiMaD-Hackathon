# Makefile
# GNU makefile for example Cahn-Hilliard model code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = $(MMSP_PATH)/include

# compilers/flags
compiler = g++
flags = -O3 -I $(incdir)

# dependencies
core = $(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp

# the program
mmsp2png: mmsp2png.cpp $(core) /usr/include/IL/devil_cpp_wrapper.hpp
	$(compiler) $(flags) -I /usr/include/IL -include il.h $< -o $@ -lz -lIL -lILU -lILUT

clean:
	rm -f mmsp2png
