#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
// #include "TDMA.h"
// #include "defineCoeffMatrix.h"
//#include "defineBoundaryCoeffMatrix.h"
#include "mesh.h"

int main()
{
    //User designed parameters:
    const double Re = 100;
    const double nu = 14.8e-6;//The kinetic viscousity of the fluid
    const double rho = 1.1691;//The density of the fluid
    const double L = 0.1;//unit length
    const double eps = 1e-7;//convergence criterion
    const int Nk = 1e+5;//iteration step
    double mu = nu*rho;
    double Um = 1.5*Re*nu/(3*L);//inlet max velocity
    double pb = 0;//back pressure
    const int N = 10;//The grid point number of unit length L
    int Nx_p = 3*(N-1)+N+6*(N-1);//grid point number along x direction for pressure
    int Ny_p = (N-1)+N+(N-1);//grid point number along y direction for pressure
    int Nx_u = Nx_p - 1;
    int Ny_u = Ny_p;
    int Nx_v = Nx_p;
    int Ny_v = Ny_p - 1;

    /*****if the solid body exists*****/
    //These are index i and j:
    double BoundaryPointLeft = 3*(N-1);
    double BoundaryPointRight = BoundaryPointLeft + (N-1);
    double BoundaryPointLow = N-1;
    double BoundaryPointUp = BoundaryPointLow + (N-1);
    /**********************************/

    /**********User defined options**********/
    const bool solidBody = true;//true: flow around a square; false: channel flow
    /****************************************/

    /*****defination of quantities*****/
    //initial condition of u and v:
    std::vector<std::vector<double>> P(Nx_p, std::vector<double> (Ny_p,0));//old field
    std::vector<std::vector<double>> Pstar(Nx_p, std::vector<double> (Ny_p,0));//new field
    std::vector<std::vector<double>> Up(Nx_p, std::vector<double> (Ny_p,0));
    std::vector<std::vector<double>> Vp(Nx_p, std::vector<double> (Ny_p,0));
    std::vector<std::vector<double>> U(Nx_u, std::vector<double> (Ny_u,0));
    std::vector<std::vector<double>> V(Nx_v, std::vector<double> (Ny_v,0));
    std::vector<std::vector<double>> gradPx(Nx_u, std::vector<double> (Ny_u, 0));
    std::vector<std::vector<double>> gradPy(Nx_v, std::vector<double> (Ny_v, 0));
    //predicted velocity u* and v*:
    std::vector<std::vector<double>> Ustar(Nx_u, std::vector<double> (Ny_u,0));
    std::vector<std::vector<double>> Vstar(Nx_v, std::vector<double> (Ny_v,0));
    /*********************************/

    /************mesh construction**********/
    mesh Mesh1(Nx_p, Ny_p, Nx_u, Ny_u, Nx_v, Ny_v, L, pb, U, V, P, Pstar);
    double dx = Mesh1.Dx();
    double dy = Mesh1.Dy();
    double _dx = 1/dx;
    double _dy = 1/dy;
    //Mesh check:
    std::cout << "\n  Mesh Information: \n" 
                << "dx = " << dx << "\n"
                << "dy = " << dy << "\n"
                << "\n  nodes number:\n" 
                << "Nx_p = " << Nx_p << "\n"
                << "Ny_p = " << Ny_p << "\n"
                << "Nx_u = " << Nx_u << "\n"
                << "Ny_u = " << Ny_u << "\n"
                << "Nx_v = " << Nx_v << "\n"
                << "Ny_v = " << Ny_v << "\n"
                <<std::endl;

    /************Including the boundary conditions**********************/
    Mesh1.defineTheBoundary(BoundaryPointLeft,
                            BoundaryPointRight, 
                            BoundaryPointLow, 
                            BoundaryPointUp,
                            L,
                            nu,
                            rho,
                            Um,
                            Ny_p,
                            Nx_u,
                            Ny_u,
                            Nx_v,
                            Ny_v,
                            P,
                            Up,
                            Vp,
                            U,
                            V,
                            !solidBody);
    
    //Inlet velocity with Re =100:
    std::cout << "**********check inlet velocity**********"<< std::endl;
    for(int j = 0; j<Ny_p; ++j)
    {
        std::cout << "Up [0]" << "[" << j << "] = " << Up[0][j] << std::endl; 
    }
    std::cout << "**********check inlet pressure**********"<< std::endl;
    for(int j = 0; j<Ny_p; ++j)
    {
        std::cout << "P [0]" << "[" << j << "] = " << P[0][j] << std::endl; 
    }

    /*******************************************************************/
    /************projection method time loop:************/
    double t=0;
    const double dt = 0.001;//(0.2*dx)/Um;//time step
    double _dt = 1.0/dt;
    double endTime = 5000;
    int writeInterval = 200;
    double x, y;
    std::cout << "**********Check the time step**********\n"
            << "dt = " << dt << std::endl;

    for(int n=0; n < endTime; n++)
    {
        t += dt;
        
        //predictor:
        #include "UEqn.h"

        //Possion solver:
        // #include "ADI.h"
        #include "PossionIteration.h"       
        //corrector:
        #include "corrector.h"

        Mesh1.defineTheBoundary(BoundaryPointLeft,
                                BoundaryPointRight, 
                                BoundaryPointLow, 
                                BoundaryPointUp,
                                L,
                                nu,
                                rho,
                                Um,
                                Ny_p,
                                Nx_u,
                                Ny_u,
                                Nx_v,
                                Ny_v,
                                P,
                                Up,
                                Vp,
                                U,
                                V,
                                !solidBody);
        std::cout << "  Time step end!  " << std::endl;
    
        //I/O:
        if(n%writeInterval==0)
        {
            std::ofstream fout;
            fout.open("./pressure/P_" + std::to_string(n) + ".tec");
            if(fout.is_open())
            {
                fout << "VARIABLES =\"x\" \n"
                        << "\"y\" \n"
                        << "\"P\" \n";

                fout << "ZONE T = \"Rank" << 1<<"\" \n"
                        << "I ="<<Nx_p<<"\n"
                        << "J ="<<Ny_p<<"\n";
                fout << "DATAPACKING = POINT\n";
                y =  -dy;
                for(int j=0; j<Ny_p; ++j)
                {   

                    x=-dx;
                    y=y+dy;
                    for(int i=0; i<Nx_p; ++i)
                    {

                        x=x+dx;
                        fout << x << " " << y << " " << P[i][j] << "\n";

                    }
                }
            }
            fout.close();
        }
        
        //I/O:
        if(n%writeInterval==0)
        {
            std::ofstream fout;
            fout.open("./velocity/U_" + std::to_string(n) + ".tec");
            if(fout.is_open())
            {
                fout << "VARIABLES =\"x\" \n"
                        << "\"y\" \n"
                        << "\"U\" \n";

                fout << "ZONE U = \"Rank" << 1<<"\" \n"
                        << "I ="<<Nx_u<<"\n"
                        << "J ="<<Ny_u<<"\n";
                fout << "DATAPACKING = POINT\n";
                y =  -dy;
                for(int j=0; j<Ny_u; ++j)
                {   

                    x=-dx;
                    y=y+dy;
                    for(int i=0; i<Nx_u; ++i)
                    {

                        x=x+dx;
                        fout << x << " " << y << " " << U[i][j] << "\n";

                    }
                }
            }
            fout.close();
        }

    }
    
    
}