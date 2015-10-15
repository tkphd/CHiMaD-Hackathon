// File:    mmsp2png.cpp
// Purpose: reads MMSP grid and writes PNG image file
// Output:  grayscale portable network graphics
// Depends: MMSP, DevIL image library, zlib

// Questions/Comments to kellet@rpi.edu (Trevor Keller)

// DevIL usage after http://bobobobo.wordpress.com/2009/03/02/how-to-use-openil-to-generate-and-save-an-image/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <zlib.h>
#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>
#include "devil_cpp_wrapper.hpp"

#include "MMSP.hpp"

void writePNG(const int w, const int h, const int bpp, unsigned char* imData, const char* filename);

//template <int dim, class U, typename T>
//void grid2png(const int field, const MMSP::grid<dim,U<T> >& grid, const char*& filename);

//template <class T>
//template <class U>
//void grid2png(const MMSP::grid<2,U<T> >& grid, const char*& filename);

int main(int argc, char* argv[])
{
  if ( argc != 3 )
  {
    std::cout << "Usage: " << argv[0] << " data.dat output.png\n";
    return ( 1 );
  }

  // file open error check
  std::ifstream input(argv[1]);
  if (!input) {
    std::cerr<<"File input error: could not open "<<argv[1]<<".\n\n";
    exit(-1);
  }

  // read data type
  std::string type;
  getline(input,type,'\n');

  // grid type error check
  if (type.substr(0,4)!="grid") {
    std::cerr<<"File input error: file does not contain grid data."<<std::endl;
    exit(-1);
  }

	// parse data type
  bool uchar_type = (type.find("unsigned char") != std::string::npos);
  bool float_type = (type.find("float") != std::string::npos);
  bool double_type = (type.find("double") != std::string::npos);
  bool long_double_type = (type.find("long double") != std::string::npos);

  bool scalar_type = (type.find("scalar") != std::string::npos);
  bool vector_type = (type.find("vector") != std::string::npos);
  bool sparse_type = (type.find("sparse") != std::string::npos);

  if (not uchar_type	 and	not float_type	and
      not double_type  and  not long_double_type) {
    std::cerr << "File input error: unknown grid data type." << std::endl;
    exit(-1);
  }

	if (not vector_type and
			not sparse_type and
			not scalar_type) {
		scalar_type=true;
	}

  // read grid dimension
 	int dim;
  input>>dim;
 	if (dim!=2) {
   	std::cerr<<"ERROR: Expected 2D input.\n"<<std::endl;
 	  exit(-1);
  }

	// read number of fields
  int fields;
  input >> fields;

  // read grid sizes
  int x0[3] = {0, 0, 0};
  int x1[3] = {0, 0, 0};
  for (int i = 0; i < dim; i++)
    input >> x0[i] >> x1[i];

  // read cell spacing
  float dx[3] = {1.0, 1.0, 1.0};
  for (int i = 0; i < dim; i++)
    input >> dx[i];

  // ignore trailing endlines
  input.ignore(10, '\n');

	int width=x1[0]-x0[0];
	int height=x1[1]-x0[1];
	unsigned int theSize=height*width;
	unsigned char* buffer = new unsigned char[theSize];

  if (scalar_type) {
		  if (uchar_type) {
	  		MMSP::grid<2,unsigned char> grid(argv[1]);
	 			for (int n=0; n<nodes(grid); ++n)
					buffer[n] = grid(n);
			}
		  else if (float_type) {
	  		MMSP::grid<2,float> grid(argv[1]);
			  float max=std::numeric_limits<float>::min();
			  float min=std::numeric_limits<float>::max();
	 			for (int n=0; n<nodes(grid); ++n) {
	 				float i=grid(n);
	 				if (i>pow(MMSP::g1(grid,0)-MMSP::g0(grid,0),dim)-1) continue;
 					else if (i>max) max=i;
			 		else if (i<min) min=i;
			 	}
			 	std::cout<<"Rescaling for φ in ["<<min<<", "<<max<<"]."<<std::endl;
		 		for (int n=0; n<nodes(grid); ++n)
					buffer[n] = (grid(n)>pow(MMSP::g1(grid,0)-MMSP::g0(grid,0),dim)-1)?255:255*(grid(n)-min)/(max-min);
			} else if (double_type) {
	  		MMSP::grid<2,double> grid(argv[1]);
			  float max=std::numeric_limits<double>::min();
			  float min=std::numeric_limits<double>::max();
	 			for (int n=0; n<nodes(grid); ++n) {
	 				double i=grid(n);
	 				if (i>pow(MMSP::g1(grid,0)-MMSP::g0(grid,0),dim)-1) continue;
 					else if (i>max) max=i;
			 		else if (i<min) min=i;
			 	}
			 	std::cout<<"Rescaling for φ in ["<<min<<", "<<max<<"]."<<std::endl;
		 		for (int n=0; n<nodes(grid); ++n)
					buffer[n] = (grid(n)>pow(MMSP::g1(grid,0)-MMSP::g0(grid,0),dim)-1)?255:255*(grid(n)-min)/(max-min);
		  } else {
		    std::cerr<<"File input error: png from "<<type<<" not implemented."<<std::endl;
    		delete [] buffer; buffer=NULL;
    		exit(-1);
    	}
  } else if (vector_type) {
  	if (fields>2) {
  		int field=0;
		  if (float_type) {
	  		MMSP::grid<2,MMSP::vector<float> > grid(argv[1]);
			  float max=std::numeric_limits<float>::min();
			  float min=std::numeric_limits<float>::max();
	 			for (int n=0; n<nodes(grid); ++n) {
 					if (grid(n)[field]>max) max=grid(n)[field];
			 		else if (grid(n)[field]<min) min=grid(n)[field];
			 	}
			 	std::cout<<"Rescaling for φ in ["<<min<<", "<<max<<"]."<<std::endl;
		 		for (int n=0; n<nodes(grid); ++n)
					buffer[n] = 255*(grid(n)[field]-min)/(max-min);
			} else if (double_type) {
	  		MMSP::grid<2,MMSP::vector<double> > grid(argv[1]);
			  float max=std::numeric_limits<double>::min();
			  float min=std::numeric_limits<double>::max();
	 			for (int n=0; n<nodes(grid); ++n) {
 					if (grid(n)[field]>max) max=grid(n)[field];
			 		else if (grid(n)[field]<min) min=grid(n)[field];
			 	}
			 	std::cout<<"Rescaling for φ in ["<<min<<", "<<max<<"]."<<std::endl;
		 		for (int n=0; n<nodes(grid); ++n)
					buffer[n] = 255*(grid(n)[field]-min)/(max-min);
		  } else {
		    std::cerr<<"File input error: png from "<<type<<" not implemented."<<std::endl;
    		delete [] buffer; buffer=NULL;
    		exit(-1);
    	}
		} else {
		  if (float_type) {
	  		MMSP::grid<2,MMSP::vector<float> > grid(argv[1]);
			 	for (int n=0; n<nodes(grid); ++n) {
			  	float sum=0;
  				for (int i=0; i<fields; ++i)
		  			sum += grid(n)[i]*grid(n)[i];
					buffer[n] = 255*sqrt(sum);
				}
 			} else if (double_type) {
	  		MMSP::grid<2,MMSP::vector<double> > grid(argv[1]);
			 	for (int n=0; n<nodes(grid); ++n) {
			  	float sum=0;
  				for (int i=0; i<fields; ++i)
		  			sum += grid(n)[i]*grid(n)[i];
					buffer[n] = 255*sqrt(sum);
				}
		  } else {
		    std::cerr<<"File input error: png from "<<type<<" not implemented."<<std::endl;
    		delete [] buffer; buffer=NULL;
    		exit(-1);
    	}
		}
  } else if (sparse_type) {
	  if (float_type) {
		 	std::string altfile(argv[2]);
		 	altfile.resize(altfile.length()-4);
		 	altfile.append("_x.png");
	  	// Image |{φ}|
	  	MMSP::grid<2,MMSP::sparse<float> > grid(argv[1]);
		  float max=std::numeric_limits<float>::min();
		  float min=std::numeric_limits<float>::max();
 			for (int n=0; n<nodes(grid); ++n) {
  			int nonzero = MMSP::length(grid(n));
		  	float sum=0;
 				for (int i=0; i<nonzero; ++i)
	  			sum += grid(n).value(i)*grid(n).value(i);
	  		if (sqrt(sum)>pow(MMSP::g1(grid,0)-MMSP::g0(grid,0),dim)-1) continue;
				else if (sqrt(sum)>max) max=sqrt(sum);
		 		else if (sqrt(sum)<min) min=sqrt(sum);
		 	}
		 	std::cout<<"Rescaling for |φ| in ["<<min<<", "<<max<<"]."<<std::endl;
		 	for (int n=0; n<nodes(grid); ++n) {
  			int nonzero = MMSP::length(grid(n));
		  	float sum=0;
  			for (int i=0; i<nonzero; ++i)
		  		sum += grid(n).value(i)*grid(n).value(i);
				buffer[n] = (sqrt(sum)>pow(MMSP::g1(grid,0)-MMSP::g0(grid,0),dim)-1)?255:255*(sqrt(sum)-min)/(max-min);
 			}
		} else if (double_type) {
	  	// Image |{φ}|
	  	MMSP::grid<2,MMSP::sparse<double> > grid(argv[1]);
		  float max=std::numeric_limits<double>::min();
		  float min=std::numeric_limits<double>::max();
 			for (int n=0; n<nodes(grid); ++n) {
  			int nonzero = MMSP::length(grid(n));
		  	float sum=0;
 				for (int i=0; i<nonzero; ++i)
	  			sum += grid(n).value(i)*grid(n).value(i);
	  		if (sqrt(sum)>pow(MMSP::g1(grid,0)-MMSP::g0(grid,0),dim)-1) continue;
				else if (sqrt(sum)>max) max=sqrt(sum);
		 		else if (sqrt(sum)<min) min=sqrt(sum);
		 	}
		 	std::cout<<"Rescaling for |φ| in ["<<min<<", "<<max<<"]."<<std::endl;
		 	for (int n=0; n<nodes(grid); ++n) {
  			int nonzero = MMSP::length(grid(n));
		  	float sum=0;
  			for (int i=0; i<nonzero; ++i)
		  		sum += grid(n).value(i)*grid(n).value(i);
				buffer[n] = (sqrt(sum)>pow(MMSP::g1(grid,0)-MMSP::g0(grid,0),dim)-1)?255:255*(sqrt(sum)-min)/(max-min);
 			}
		}
  } else {
    std::cerr<<"File input error: png from "<<type<<" not implemented."<<std::endl;
    delete [] buffer; buffer=NULL;
    exit(-1);
  }
	writePNG(width, height, 1, buffer, argv[2]);
 	delete [] buffer; buffer=NULL;

  return 0;
}

void writePNG(const int w, const int h, const int bpp, unsigned char* imData, const char* filename)
{
  // Initialize image
 	ilInit();
 	ILenum Error;
  ILuint imageID = ilGenImage() ;
 	ilBindImage(imageID);
  ilTexImage(w, h, 1, bpp, IL_LUMINANCE, IL_UNSIGNED_BYTE, imData);
 	Error = ilGetError();
  if (Error!=IL_NO_ERROR) std::cout<<"Error making image: "<<iluErrorString(Error)<<std::endl;
 	ilEnable(IL_FILE_OVERWRITE);
  ilSave( IL_PNG, filename) ;
 	Error = ilGetError();
  if (Error!=IL_NO_ERROR) std::cout<<"Error saving image: "<<iluErrorString(Error)<<std::endl;
}
