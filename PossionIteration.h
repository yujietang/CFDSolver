/*
Solving the Possion equation;
The Iteration Method of the Gauss-Seidel Algorithm is used; 
k: the iteration step;
*/
for(int j=0; j<Ny_p; ++j)
{
    for(int i=0; i<Nx_p; ++i)
    {
        std::cout << "/******************Scan Point*******************/\n" << std::endl;
        std::cout <<"point i = "<< i << std::endl;
        std::cout <<"point j = "<< j <<"\n"<< std::endl;

        if(//inlet**********
            i==0 && j>0 && j<(Ny_p-1)
        )
        {  
            Pstar[i][j] = P[i+1][j] 
                        - (dx*mu)*((_dx*_dx)*(0.5*(Ustar[i+1][j]+Ustar[i][j])
                        - 2*Ustar[i][j] + Up[i][j])
                        + (_dy*_dy)*(Vstar[i][j]-2*Vp[i][j]+Vstar[i][j-1]));
            P[i][j] = Pstar[i][j];     
        }
        else if(//outlet**********       
            (i==Nx_p-1) && j>0 && j<(Ny_p-1)
        )
        {

            Pstar[i][j] = pb;
            P[i][j] = Pstar[i][j];
        }
        else if(//upper wall**********
            j == Ny_p-1 && i>=0 && i<Nx_p
        )
        {
            
        for(int k=0; k<Nk; ++k)
            {
                Pstar[i][j] = P[i][j-1] 
                        + (dy*mu*_dy*_dy)*(Vp[i][j]-2*Vstar[i][j-1]
                        + 0.5*(Vstar[i][j-1]+Vstar[i][j-2]));
                if((abs(1-(Pstar[i][j]/P[i][j]))) < eps)
                {
                    break;
                }
                P[i][j] = Pstar[i][j];
            }
        }
        else if(//lower wall**********
            j==0 && i>=0 && i<Nx_p
        )
        {
            
            for(int k=0; k<Nk; ++k)
            {
                Pstar[i][j] = P[i][j+1] 
                        - (dy*mu*_dy*_dy)*(Vp[i][j]-2*Vstar[i][j]
                        + 0.5*(Vstar[i][j]+Vstar[i][j+1]));
                if((abs(1-(Pstar[i][j]/P[i][j]))) < eps)
                {
                    break;
                }
                P[i][j] = Pstar[i][j];
            }
        }
        else
        {
            for(int k=0; k<Nk; ++k)
            {
                std::cout << "Internal Field Iteration: k = " << k << std::endl;
                Pstar[i][j] = (1/(2*_dx*_dx+2*_dy*_dy))
                        *((1*_dx*_dx)*P[i-1][j] + (1*_dy*_dy)*P[i][j-1] 
                        +(1*_dx*_dx)*P[i+1][j] + (1*_dy*_dy)*P[i][j+1]
                        -(rho*_dt*_dx)*(Ustar[i][j]-Ustar[i-1][j])
                        -(rho*_dt*_dy)*(Vstar[i][j]-Vstar[i][j-1]));
                if((abs(1-(Pstar[i][j]/P[i][j]))) < eps)
                {
                    break;
                }
                P[i][j] = Pstar[i][j];
            }
        }
        //Debug
        std::cout << "Pressure Field:\n"
                    << "P [" << i << "," << j << "]= " << P[i][j] << std::endl; 
        
    }
}
//4 Corner points:
P[0][0] = 0.5*(P[1][0]+P[0][1]);
P[0][Ny_p-1] = 0.5*(P[0][Ny_p-2]+P[1][Ny_p-1]);
P[Nx_p-1][0] = 0.5*(P[Nx_p-2][0]+P[Nx_p-1][1]);
P[Nx_p-1][Ny_p-1] = 0.5*(P[Nx_p-1][Ny_p-2]+P[Nx_p-2][Ny_p-1]);
/**********Substitude the Pressure Gradient**********/
//For points storing U:
for(int i=0; i<Nx_u; ++i)
{
    for(int j=0; j<Ny_u; ++j)
    {
        gradPx[i][j] = _dx*(P[i+1][j]-P[i][j]);
    }
}

//For points storing V:
for(int i=0; i<Nx_v; ++i)
{
    for(int j=0; j<Ny_v; ++j)
    {
        gradPy[i][j] = _dx*(P[i][j+1]-P[i][j]);
    }
}

