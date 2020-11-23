/***********************************************
The mesh is a Cartesian staggered grid!
The pressure is stored at the center of the grid element;
The velocity vecotr is stored at the sides of the the grid element;
***********************************************/
class mesh
{
    public:
        //Mesh construction:
        mesh(const int Nx_p, 
            const int Ny_p,
            const int Nx_u,
            const int Ny_u,
            const int Nx_v,
            const int Ny_v, 
            const double L,
            const double pb,
            std::vector<std::vector<double>>& U,
            std::vector<std::vector<double>>& V,
            std::vector<std::vector<double>>& P,
            std::vector<std::vector<double>>& Pstar)
        {
        /*****initialize the discretized parameters*****/
            dx = (10*L)/(Nx_p-1);
            dy = (3*L)/(Ny_p-1);

        /*****initial conditions of internal field*****/
            for(int i=0; i<Nx_u; ++i)
            {
                for(int j=1; j<(Ny_u-1); ++j)
                {
                    U[i][j] = 0;
                }
            }
            for(int i=1; i<(Nx_v-1); ++i)
            {
                for(int j=0; j<Ny_v; ++j)
                {
                    V[i][j] = 0;
                }
            }
            for(int i=0; i<Nx_p; ++i)
            {
                for(int j=0; j<Ny_p; ++j)
                {
                    P[i][j] = pb;
                    Pstar[i][j] = P[i][j];
                    // std::cout << "Debug: P[" << i << "]" << "[" << j << "] = " << P[i][j] << std::endl;
                }
            }

            // std::cout << "dx = " << dx << std::endl;
            // std::cout << "dy = " << dy << std::endl;
        }
        /**********Boundary Identification**********/
        void defineTheBoundary(
            const int BoundaryPointLeft, 
            const int BoundaryPointRight, 
            const int BoundaryPointLow, 
            const int BoundaryPointUp,
            const double L,
            const double nu,
            const double rho,
            const double Um,
            const int Ny_p,
            const int Nx_u,
            const int Ny_u,
            const int Nx_v,
            const int Ny_v,
            std::vector<std::vector<double>>& P,
            std::vector<std::vector<double>>& Up,
            std::vector<std::vector<double>>& Vp,
            std::vector<std::vector<double>>& U,
            std::vector<std::vector<double>>& V,
            const bool solidBodyExist)
        {
            /**********Initial Condition**********/
            //The velocity has a parabolic profile
            double yp = 0;//y position in Cartesian coordinate
            double H = (Ny_p-1)*dy;
            std::cout << "********** Umax **********" << Um << std::endl;

            for(int jj=1; jj<Ny_p-1; ++jj){
                yp = jj*dy;
                Up[0][jj] = 2*Um*(1-yp/(3*L))*yp;//yp*yp*(-Um/(2.25*L*L))+yp*(-3*L*(-Um/(2.25*L*L)));
                Vp[0][jj] = 0;
                P[0][jj] = P[1][jj] + (Up[0][jj]*2*nu*rho)/(yp*(H-yp));
                std::cout << "P [0," << jj <<"]= " << P[0][jj] << std::endl; 
            }
            /*************************************/
            if(solidBodyExist)
            {
                //The solid square inside the flow channel case:
                //internal region of solid body:
                for(int i=BoundaryPointLeft; i<=(BoundaryPointRight - 1); ++i)
                {
                    for(int j=BoundaryPointLow; j<=BoundaryPointUp; ++j)
                    {
                        U[i][j] = 0;
                    }
                }
                for(int i=BoundaryPointLeft; i<=BoundaryPointRight; ++i)
                {
                    for(int j=BoundaryPointLow; j<=(BoundaryPointUp-1); ++j)
                    {
                        V[i][j] = 0;
                    }
                }
                /***************Velocity Boundary Condition***************/
                //left and right wall of solid square:
                for(int j=BoundaryPointLow; j<=BoundaryPointUp; ++j)
                {
                    //left surface:
                    Up[BoundaryPointLeft][j] = 0;
                    //right surface:
                    Up[BoundaryPointRight][j] = 0;
                    /*********************************************************/
                }
                for(int j=BoundaryPointLow; j<=(BoundaryPointUp-1); ++j)
                {
                    //Left surface:
                    V[BoundaryPointLeft][j] = 0;
                    //Right surface:
                    V[BoundaryPointRight][j] = 0;
                }

                //upper and lower wall of solid square:
                for(int i=BoundaryPointLeft; i<=(BoundaryPointRight-1); ++i)
                {
                    /***************Velocity Boundary Condition***************/
                    //Upper surface:
                    U[i][BoundaryPointUp] = 0;
                    //Lower surface:
                    U[i][BoundaryPointLow] = 0;
                    /*********************************************************/
                }
                for(int i=BoundaryPointLeft; i<BoundaryPointRight; ++i)
                {
                    //Upper surface:
                    Vp[i][BoundaryPointUp] = 0;
                    //Lower surface:
                    Vp[i][BoundaryPointLow] = 0;
                }
            }
            else
            {
                for(int i=0; i<Nx_u; ++i)
                {
                    U[i][0] = 0;
                    U[i][Ny_u-1] = 0;
                }
                for(int j=0; j<Ny_v; ++j)
                {
                    V[0][j] = 0;
                    V[Nx_v-1][j] = 0;
                }
            }
        }
        
        //Member functions:
        double Dx() const
        {
            return dx;
        }
        double Dy() const
        {
            return dy;
        }
    private:
        /**********Grid discretized size**********/
        double dx;//cell size delta x
        double dy;//cell size delta y
};