// cahn-hilliard.hpp
// Example code for 2D Cahn-Hilliard model problem
// Questions/comments to trevor.keller@nist.gov (Trevor Keller)

#ifndef CAHNHILLIARD_UPDATE
#define CAHNHILLIARD_UPDATE
#include"MMSP.hpp"          // pulls in the MMSP back-end (grid, data classes, etc.)
#include<cmath>
#include"cahn-hilliard.hpp" // defines the grid dimension (2) and type (double)

// Hard-code some system parameters and numerical constants at global scope for convenience

const double deltaX = 1.0;
const double Ca = 0.05;
const double Cb = 0.95;
const double Cm = 0.5*(Ca + Cb);
const double A = 2.0;
const double B = A/((Ca-Cm)*(Ca-Cm));
const double D = 2.0/(Cb-Ca);
const double K = 2.0;
const double CFL = 0.25;
const double dt = std::pow(deltaX, 4)*CFL/(32.0*D*K);

double energydensity(const double& C)
{
    return -0.5*A*pow(C-Cm,2) + 0.25*B*pow(C-Cm,4) + 0.25*Ca*pow(C-Ca,4) + 0.25*Cb*pow(C-Cb,4);
}

double dfdc(const double& C)
{
    return -A*(C-Cm) + B*pow(C-Cm, 3) + Ca*pow(C-Ca, 3) + Cb*pow(C-Cb, 3);
}


namespace MMSP {

void generate(int dim, const char* filename)
{
    if (dim!=2) {
        std::cerr<<"ERROR: Initial condition is 2-D, only!"<<std::endl;
        std::exit(-1);
    }

    int rank=0;
    #ifdef MPI_VERSION
    rank = MPI::COMM_WORLD.Get_rank();
    #endif                                                                                                                                               

    const double q[2] = {0.1*std::sqrt(2.0), 0.1*std::sqrt(3.0)};                                                                                            
                                                                                                                                                             
    if (dim==2) {                                                                                                                                        
        MMSP::grid<2,double> grid(1,0,200,0,200);                                                                                                    
        for (int d=0; d<dim; d++){                                                                                                                   
            dx(grid,d) = deltaX;                                                                                                                 
            if (MMSP::x0(grid,d)==MMSP::g0(grid,d))                                                                                              
                MMSP::b0(grid,d) = Neumann; // enumerated in MMSP.utility.hpp                                                                    
            else if (MMSP::x1(grid,d)==MMSP::g1(grid,d))                                                                                         
                MMSP::b1(grid,d) = Neumann; // enumerated in MMSP.utility.hpp                                                                    
        }                                                                                                                                            

        for (int i=0; i<nodes(grid); i++) {                                                                                                          
            MMSP::vector<int> x = position(grid,i);                                                                                              
            grid(x) = 0.45 + 0.01 * std::cos(x[0]*dx(grid,0)*q[0] + x[1]*dx(grid,1)*q[1]);                                                       
        }                                                                                                                                            

        #ifdef MPI_VERSION                                                                                                                           
        MPI::COMM_WORLD.Barrier();                                                                                                                   
        #endif                                                                                                                                       
        output(grid,filename);                                                                                                                       
        if (rank==0)                                                                                                                                 
            std::cout<<"Timestep is "<<dt<<" (Co="<<CFL<<')'<<std::endl;                                                                         
    }                                                                                                                                                    
}                                                                                                                                                            
                                                                                                                                                             
template <int dim, typename T>                                                                                                                               
void update(MMSP::grid<dim,T>& grid, int steps)                                                                                                              
{                                                                                                                                                            
    int rank=0;                                                                                                                                          
    #ifdef MPI_VERSION                                                                                                                                   
    rank = MPI::COMM_WORLD.Get_rank();                                                                                                                   
    #endif                                                                                                                                               

    // Make sure the grid spacing is correct                                                                                                             
    for (int d=0; d<dim; d++) {                                                                                                                          
        dx(grid,d) = deltaX;                                                                                                                         
        if (MMSP::x0(grid,d)==MMSP::g0(grid,d))                                                                                                      
            MMSP::b0(grid,d) = Neumann; // enumerated in MMSP.utility.hpp                                                                            
        else if (MMSP::x1(grid,d)==MMSP::g1(grid,d))                                                                                                 
            MMSP::b1(grid,d) = Neumann; // enumerated in MMSP.utility.hpp                                                                            
    }                                                                                                                                                        
                                                                                                                                                             
    ghostswap(grid);                                                                                                                                     
                                                                                                                                                             
    // Let us be absolutely explicit about BCs here.                                                                                                          
    MMSP::grid<dim,T> update(grid);
    for (int d=0; d<dim; d++) {
        dx(update,d) = deltaX;
        if (MMSP::x0(update,d)==MMSP::g0(update,d))
            MMSP::b0(update,d) = Neumann; // enumerated in MMSP.utility.hpp
        else if (MMSP::x1(update,d)==MMSP::g1(update,d))
            MMSP::b1(update,d) = Neumann; // enumerated in MMSP.utility.hpp
    }

    MMSP::grid<dim,T> temp(grid);
    for (int d=0; d<dim; d++) {
        dx(temp,d) = deltaX;
        if (MMSP::x0(temp,d)==MMSP::g0(temp,d))
            MMSP::b0(temp,d) = Neumann; // enumerated in MMSP.utility.hpp
        else if (MMSP::x1(temp,d)==MMSP::g1(temp,d))
            MMSP::b1(temp,d) = Neumann; // enumerated in MMSP.utility.hpp
    }


    for (int step=0; step<steps; step++) {
        for (int i=0; i<nodes(grid); i++) {
            MMSP::vector<int> x = position(grid,i);
            double c = grid(x);
            temp(i) = dfdc(c) - K*laplacian(grid,x);
        }

        ghostswap(temp);

        double energy = 0.0;
        double mass = 0.0;
        for (int i=0; i<nodes(grid); i++) {
                MMSP::vector<int> x = position(grid,i);
                update(x) = grid(x)+dt*D*laplacian(temp,x);
                energy += dx(grid)*dy(grid)*energydensity(update(x));
                mass += dx(grid)*dy(grid)*update(x);
        }
        #ifdef MPI_VERSION
        MPI::COMM_WORLD.Barrier();
        double myEnergy = energy;
        double myMass = mass;
        MPI::COMM_WORLD.Allreduce(&myEnergy, &energy, 1, MPI_DOUBLE, MPI_SUM);
        MPI::COMM_WORLD.Allreduce(&myMass, &mass, 1, MPI_DOUBLE, MPI_SUM);
        #endif
        #ifndef DEBUG
        if (rank==0)
            std::cout<<energy<<'\t'<<mass<<'\n';
        #endif

        swap(grid,update);
        ghostswap(grid);
    }
    #ifndef DEBUG
    if (rank==0)
        std::cout<<std::flush;
    #endif
}

} // MMSP
#endif

#include"MMSP.main.hpp"