//
//  main.cpp
//
//  Created by Emad Ghalenoei in 2020.
//  Copyright (c) 2020 Emad.Ghalenoei. All rights reserved.

#include <fstream> // for file access
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <sstream>
#include <string>
#include <cmath> // For Absolute function
//#include "matrix.h" // Matrix Class
#include <random>
#include <chrono>
#include "matrix.hpp" // Matrix Class
#include <tuple>


#define PI 3.14159265

using namespace std;

int main()
{

srand(time(NULL));


    Matrix kernel=Matrix("kernel.txt"); // you need to define the kernel or load it separately
    //Matrix CM=Matrix("CM.txt");
    Matrix dg_load=Matrix("data.txt");
    Matrix C_original=Matrix("C_original.txt");
    int num_obs = dg_load.getRows();

    Matrix x_s     = dg_load.MatExtraction(1, num_obs, 1, 1);
    Matrix dg_obs  = dg_load.MatExtraction(1, num_obs, 2, 2);
    Matrix dg_true = dg_load.MatExtraction(1, num_obs, 3, 3);
    const double true_e=dg_load(0,3);
    // double density = -0.4;   // Unit(kg/m^3)
    const double density = -0.4 * 1000;  // Unit(kg/m^3)

    Matrix x=x_s;
    const double Z_model_space = 100;      // difference betweeen each depth layer (m)
    const double Z1=0;                     // first depth at model space (m)
    const double Z2=10000;                 // final depth at model space (m)
    //int num_layer=100;                 // number of z layer
    //Matrix z (num_layer,1);
    //z = Matrix::linespace(z(0,0), Z2, num_layer);
    Matrix z = Matrix::colon(Z1+(Z_model_space/2),Z_model_space,Z2);
    const int NX = x.getRows();
    const int NZ = z.getRows();
    Matrix X(NZ,NX);
    Matrix Z(NZ,NX);
    // tie(X, Z) = Matrix::meshgrid(x, z, X, Z,NX,NZ);
    Matrix::meshgrid(x, z, X, Z,NX,NZ);

    Matrix true_noise=dg_obs-dg_true;
    double True_LOGL = Matrix::Log_Likelihood_full_C(true_noise,C_original);

    // initializing
    int n=6;
    const double x_min=X.Min();
    const double x_max=X.Max();
    const double z_min=Z.Min();
    const double z_max=Z.Max();
    Matrix xr =  Matrix::RandomDouble(n,1)*(x_max-x_min) + x_min;
    Matrix zr =  Matrix::RandomDouble(n,1)*(z_max-z_min) + z_min;
    Matrix Bc = Matrix::join_in_right({xr, zr, Matrix::RandomDouble(n,1).rounding()});
    Matrix cc = Matrix::XZ_model(Bc, X, Z,NX,NZ);
    Matrix dc = Matrix::FM_fast(cc,kernel,density);
    Matrix res = dg_obs-dc;
    double ec = (Matrix::rand01() * (0.1-0.01)) + 0.01 ;
    Matrix dg_obs_abs = dg_obs.absolute();
    //Matrix C = Matrix::e2C(ec,dg_obs_abs,num_obs);
    double LogLc = Matrix::Log_Likelihood_initial(res,dg_obs_abs,ec);
    int Nmin=4;
    int Nmax=10;
    Matrix F_TAB (1,3);
    F_TAB (0,0)=LogLc ;
    F_TAB (0,1)=Bc.getRows() ;
    F_TAB (0,2)=ec ;
    Matrix zeros(1,(Nmax-Bc.getRows())*3);
    Matrix tab = Matrix::join_in_right({F_TAB,Bc.vectorize().transpose(),zeros});
    Matrix T_SA = Matrix::logspace(0.7,0,15);
    int SA_ITE=1000;
    //SA_TAB=Matrix::MCMC_0(X,Z,dg_obs,kernel,SA_TAB,density,Nmin,Nmax,T_SA(i,0),k,True_LOGL);


    for (unsigned i=0; i<T_SA.getRows() ; ++i)
    {
        for (int k=1; k<=SA_ITE; ++k)
        {
            Matrix::MCMC_0(X,Z,x_min,x_max,z_min,z_max,NX,NZ,dg_obs,dg_obs_abs,num_obs,kernel,tab,density,Nmin,Nmax,T_SA(i,0),k,True_LOGL);

        }
    }



    return 0;
}
