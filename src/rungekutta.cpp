/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: rungekutta.cpp
// Purpose: implementation of various Runge-Kutta methods

#ifndef RUNGEKUTTA_CPP
#define RUNGEKUTTA_CPP

#include "rungekutta.hpp"

#include <iostream> // cout endl;
//#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE

#define ONESIXTH 0.16666666666666666666666666666666666666666666666666666666

double TEM_NS::marcusRK4( 
                     double (*integrand)(const double& tt, const double& yy,
                                          const double& poissonJump), 
                     const double& dt, 
                     const double& y0, 
                     const double& t0, 
                     const double& jump)
{
   double k1, k2, k3, k4; 
   k1 = dt* integrand( t0, y0, jump);
   k2 = dt* integrand( t0 + 0.5*dt, y0 + 0.5*k1, jump);
   k3 = dt* integrand( t0 + 0.5*dt, y0 + 0.5*k2, jump);
   k4 = dt* integrand( t0 + dt, y0 + k3, jump);

   //std::cout << "t0: " << t0 << ", y0:" << y0  << ", dt: " << dt << ", jump: " << jump << ", k's: " << k1 << ", " << k2 << ", " << k3 << ", " << k4 << std::endl;

   return y0 + ONESIXTH *(k1 + 2.0*k2 + 2.0*k3 + k4 );
}

#endif
