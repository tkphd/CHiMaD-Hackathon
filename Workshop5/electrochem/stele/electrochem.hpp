// electrochem.cpp
// Benchmark program for the electrochem model using MMSP
// Questions/comments to trevor.keller@nist.gov (Trevor Keller)

std::string PROGRAM = "electrochem";
std::string MESSAGE = "electrochem model using MMSP";

typedef MMSP::grid<1,MMSP::vector<double> > GRID1D;
typedef MMSP::grid<2,MMSP::vector<double> > GRID2D;
typedef MMSP::grid<3,MMSP::vector<double> > GRID3D;

