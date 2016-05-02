// File:    linescan.cpp
// Purpose: reads MMSP grid and writes linescan down the center to CSV

// Questions/Comments to trevor.keller@gmail.com (Trevor Keller)


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <zlib.h>

#include "MMSP.hpp"
#include "energy.hpp"

int main(int argc, char* argv[])
{
  if ( argc != 3 )
  {
    std::cout << "Usage: " << argv[0] << " data.dat output.csv\n";
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

    std::ofstream of(argv[2]);
    if (!of) {
        std::cerr<<"Error: unable to open "<<argv[1]<<". Check permissions."<<std::endl;
        std::exit(-1);
    }

    if (scalar_type) {
        if (uchar_type) {
	        MMSP::grid<2,unsigned char> grid(argv[1]);
	        of<<"y,p(y)\n";
	        MMSP::vector<int> x(2, (MMSP::x1(grid)-MMSP::x0(grid))/2);
	        for (x[1]=MMSP::y0(grid); x[1]<MMSP::y1(grid); x[1]++)
	            of<<x[1]<<','<<grid(x)<<'\n';
        } else if (float_type) {
  	    	MMSP::grid<2,float> grid(argv[1]);
	        of<<"y,p(y)\n";
	        MMSP::vector<int> x(2, (MMSP::x1(grid)-MMSP::x0(grid))/2);
	        for (x[1]=MMSP::y0(grid); x[1]<MMSP::y1(grid); x[1]++)
	            of<<x[1]<<','<<grid(x)<<'\n';
    	} else if (double_type) {
  	    	MMSP::grid<2,double> grid(argv[1]);
	        of<<"y,p(y)\n";
	        MMSP::vector<int> x(2, (MMSP::x1(grid)-MMSP::x0(grid))/2);
	        for (x[1]=MMSP::y0(grid); x[1]<MMSP::y1(grid); x[1]++)
	            of<<x[1]<<','<<grid(x)<<'\n';
    	} else {
	        std::cerr<<"File input error: png from "<<type<<" not implemented."<<std::endl;
       		exit(-1);
   	    }
    } else if (vector_type) {
		if (float_type) {
	  	    MMSP::grid<2,MMSP::vector<float> > grid(argv[1]);
	        of<<"y,f";
            for (int i=0; i<MMSP::fields(grid); i++)
                of<<",p"<<i;
            of<<'\n';
	        MMSP::vector<int> x(2, (MMSP::x1(grid)-MMSP::x0(grid))/2);
	        for (x[1]=MMSP::y0(grid); x[1]<MMSP::y1(grid); x[1]++) {
	            of<<x[1]<<','<<energydensity(grid(x));
	            for (int i=0; i<MMSP::fields(grid); i++)
	                of<<','<<grid(x)[i];
	            of<<'\n';
	         }
	    } else if (double_type) {
	  	    MMSP::grid<2,MMSP::vector<double> > grid(argv[1]);
	        of<<"y,f";
            for (int i=0; i<MMSP::fields(grid); i++)
                of<<",p"<<i;
            of<<'\n';
	        MMSP::vector<int> x(2, (MMSP::x1(grid)-MMSP::x0(grid))/2);
	        for (x[1]=MMSP::y0(grid); x[1]<MMSP::y1(grid); x[1]++) {
	            of<<x[1]<<','<<energydensity(grid(x));
	            for (int i=0; i<MMSP::fields(grid); i++)
	                of<<','<<grid(x)[i];
	            of<<'\n';
	         }
        } else {
		    std::cerr<<"File input error: png from "<<type<<" not implemented."<<std::endl;
    		exit(-1);
    	}
    } else if (sparse_type) {
        if (float_type) {
	      	MMSP::grid<2,MMSP::sparse<float> > grid(argv[1]);
	        of<<"y,|p|,p1,...\n";
	        MMSP::vector<int> x(2, (MMSP::x1(grid)-MMSP::x0(grid))/2);
	        for (x[1]=MMSP::y0(grid); x[1]<MMSP::y1(grid); x[1]++) {
	            of<<x[1]<<','<<grid(x).getMagPhi();
	            for (int i=0; i<MMSP::length(grid(x)); i++)
	                of<<','<<grid(x)[i];
	            of<<'\n';
	         }
	    } else if (double_type) {
    	  	MMSP::grid<2,MMSP::sparse<double> > grid(argv[1]);
	        of<<"y,|p|,p1,...\n";
	        MMSP::vector<int> x(2, (MMSP::x1(grid)-MMSP::x0(grid))/2);
	        for (x[1]=MMSP::y0(grid); x[1]<MMSP::y1(grid); x[1]++) {
	            of<<x[1]<<','<<grid(x).getMagPhi();
	            for (int i=0; i<MMSP::length(grid(x)); i++)
	                of<<','<<grid(x)[i];
	            of<<'\n';
	         }
	}
    of.close();
  } else {
    std::cerr<<"File input error: png from "<<type<<" not implemented."<<std::endl;
    exit(-1);
  }

  return 0;
}

