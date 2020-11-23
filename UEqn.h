//Predictor to solve the U* equation:
//For one time step:

/*****dealing with the near boundary points*****/
//For U:
for(int j = 1; j < (Ny_u-1); ++j)
{
//inlet region:
    Ustar[0][j] = U[0][j] 
                + dt*(-_dx*U[0][j]*(0.5*(U[1][j]+U[0][j])-Up[0][j]) 
                - (0.5*_dy)*0.25*(V[0][j]+V[1][j]+V[0][j-1]+V[1][j-1])*(U[0][j+1]-U[0][j-1])
                + (nu*_dx*_dx)*(0.5*(U[1][j]+U[0][j])-2*U[0][j]+Up[0][j])
                + (nu*_dy*_dy)*(U[0][j+1]-2*U[0][j]+U[0][j-1]));
    // std::cout << "Ustar [0][" << j << "] = " << Ustar[0][j] << std::endl;

//outlet region:
    Ustar[Nx_u-1][j] = U[Nx_u-1][j] 
                    + dt*(-_dx*U[Nx_u-1][j]*(Up[Nx_p-1][j]-0.5*(U[Nx_u-1][j]+U[Nx_u-2][j])) 
                    - (0.5*_dy)*0.25*(V[Nx_v-1][j]+V[Nx_v-1][j-1]+V[Nx_v-2][j]+V[Nx_v-2][j-1])*(U[Nx_u-1][j+1]-U[Nx_u-1][j-1])
                    + (nu*_dx*_dx)*(0.5*(U[Nx_u-1][j]+U[Nx_u-2][j])-2*U[Nx_u-1][j]+Up[Nx_p-1][j])
                    + (nu*_dy*_dy)*(U[Nx_u-1][j+1]-2*U[Nx_u-1][j]+U[Nx_u-1][j-1]));
    // std::cout << "Ustar [Nx_u-1][" << j << "] = " << Ustar[Nx_u-1][j] << std::endl;
}

//For V:
for(int i = 1; i < (Nx_v-1); ++i)
{
//upper wall:
    Vstar[i][Ny_v-1] = V[i][Ny_v-1] 
                    + dt*(-(0.5*_dx)*0.25*(U[i][Ny_u-1]+U[i][Ny_u-2]+U[i-1][Ny_u-1]+U[i-1][Ny_u-2])*(V[i+1][Ny_v-1]-V[i-1][Ny_v-1])
                    - _dy*V[i][Ny_v-1]*(Vp[i][Ny_p-1]-0.5*(V[i][Ny_v-1]+V[i][Ny_v-2]))
                    + (nu*_dx*_dx)*(V[i+1][Ny_v-1]-2*V[i][Ny_v-1]+V[i-1][Ny_v-1])
                    + (nu*_dy*_dy)*(Vp[i][Ny_p-1]-2*V[i][Ny_v-1]+0.5*(V[i][Ny_v-1]+V[i][Ny_v-2])));
//Lower wall:
    Vstar[i][0] = V[i][0] 
                + dt*(-(0.5*_dx)*0.25*(U[i][0]+U[i][1]+U[i-1][0]+U[i-1][1])*(V[i+1][0]-V[i-1][0])
                - _dy*V[i][0]*(0.5*(V[i][0]+V[i][1])-Vp[i][0])
                + (nu*_dx*_dx)*(V[i+1][0]-2*V[i][0]+V[i-1][0])
                + (nu*_dy*_dy)*(Vp[i][0]-2*V[i][0]+0.5*(V[i][0]+V[i][1])));
    // std::cout << "Vstar [" << i << "][0] = " << Vstar[i][0] << std::endl;

}

/**********Internal Field**********/
//for U:
for(int i = 1; i < (Nx_u-1); ++i)
{
    for(int j = 1; j < (Ny_u-1); ++j)
    {
        Ustar[i][j] = U[i][j] + dt*(-(0.5*_dx)*U[i][j]*(U[i+1][j]-U[i-1][j]) 
                            - (0.5*_dy)*0.25*(V[i][j]+V[i+1][j]+V[i][j-1]+V[i+1][j-1])*(U[i][j+1]-U[i][j-1])
                            + (nu*_dx*_dx)*(U[i+1][j]-2*U[i][j]+U[i-1][j])
                            + (nu*_dy*_dy)*(U[i][j+1]-2*U[i][j]+U[i][j-1]));
        // std::cout << "Ustar [" << i << "][" << j << "]" << Ustar[i][j] << std::endl;

    }
}
//for V:
for(int i = 1; i < (Nx_v-1); ++i)
{
    for(int j = 1; j < (Ny_v-1); ++j)
    {
        Vstar[i][j] = V[i][j] + dt*(-(0.5*_dx)*0.25*(U[i][j]+U[i][j+1]+U[i-1][j]+U[i-1][j+1])*(V[i+1][j]-V[i-1][j])
                            - (0.5*_dy)*V[i][j]*(V[i][j+1]-V[i][j-1])
                            + (nu*_dx*_dx)*(V[i+1][j]-2*V[i][j]+V[i-1][j])
                            + (nu*_dy*_dy)*(V[i][j+1]-2*V[i][j]+V[i][j-1]));
        // std::cout << "Vstar [" << i << "][" << j << "]" << Vstar[i][j] << std::endl;
    }
}


