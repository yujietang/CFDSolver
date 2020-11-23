/**********internal field**********/
for(int i=0; i<Nx_u; ++i)
{
    for(int j=0; j<(Ny_u-1); ++j)
    {
        U[i][j] = Ustar[i][j] - (dt/rho)*gradPx[i][j];
    }
}

for(int i=0; i<Nx_v; ++i)
{
    for(int j=0; j<Ny_v; ++j)
    {
        V[i][j] = Vstar[i][j] - (dt/rho)*(gradPy[i][j]);
    }
}

/**********Boundary conditions(Correction)**********/
/*****Inlet*****/



