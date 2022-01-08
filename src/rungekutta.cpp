/* ----------------------------------------------------------------------
    SPF - Stochastic Phase Field
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the README file in the top-level SPF directory.
---------------------------------------------------------------------- */
// File: rungekutta.cpp

#ifndef RUNGEKUTTA_CPP
#define RUNGEKUTTA_CPP

#include "rungekutta.hpp"

#include <iostream> // cout endl;
//#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE

#ifndef ONESIXTH 
#define ONESIXTH 0.16666666666666666666666666666666666666666666666666666667
#endif

int SPF_NS::marcusIntegral( double (*marcusIntegrand)(const double& tt,
                                          const double& yy,
                                          const double& dP), 
                     const double& rk_dt,    // increment of unit interval
                     const double& y0,    // path.back()
                     const double& t0,    // t_n, t_{n+1}=t_n + dt
                     const double& dP,// Y_{k} size of jump before marcus
                     double& jumpDestination)
{
   jumpDestination = y0;
   if (( rk_dt >= 1.0) || (rk_dt <= 0))
   {
      std::cout << "Error: runge-kutta increment rk_dt >=1 or <=0" 
               << std::endl;
      return EXIT_FAILURE;
   }
   for( double tt=0.0; tt < 1.0; tt += rk_dt ) 
   {
      jumpDestination = marcusRK4( marcusIntegrand,
                                    rk_dt,
                                    jumpDestination,
                                    tt,
                                    dP);
   }
   return EXIT_SUCCESS;
}

double SPF_NS::marcusRK4( 
                     double (*integrand)(const double& tt, 
                                          const double& yy,
                                          const double& poissonJump),// Y_k
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
