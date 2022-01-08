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
// File: rand1D.cpp

#ifndef RAND1D_CPP
#define RAND1D_CPP

#include "rand1D.hpp"

#include <iostream>
using std::cout;
using std::endl;

double SPF_NS::random::exponentialWaitTime( const double& rate, const double& totalRate )
{
   double r1;
   r1 = uniform_scale( generator );
   return (1.0/totalRate)*log(1.0/r1);
}

bool SPF_NS::random::poisson_trial()
{
   if ( uniform_scale(generator) <= subinterval )
      return true;
   else
      return false;
}

int SPF_NS::random::gaussiandisplacement1D( double& xx, const double& scale)
{
   xx = scale * (2.0*uniform_scale( generator ) - 1.0);
   return EXIT_SUCCESS;
}

int SPF_NS::random::unitdisplacement1D( double& xx, const double& fwd, 
                                    const double& bkwd)
{
   if (fwd + bkwd > 1.0) 
   {
      cout << "error: displacement1D, probability >= 1.0)" << endl;
      return EXIT_FAILURE;
   }
   if ((fwd < 0.0) || (bkwd < 0.0)) 
   {
      cout << "error: displacement1D, probability < 0.0" << endl;
      return EXIT_FAILURE;
   }

   xx = uniform_scale( generator );

   if ( xx < fwd ) 
   {
      xx = 1.0;
      return EXIT_SUCCESS;
   }
   else if ( xx < (fwd + bkwd) ) 
   {
      xx = -1.0;
      return EXIT_SUCCESS;
   }
   else
   {
      xx = 0.0;
      return EXIT_SUCCESS;
   }

   return EXIT_SUCCESS;
}

bool SPF_NS::random::bernoulli( const double& pp )
{
   std::bernoulli_distribution bb( pp );
   return bb( generator );
}

//int SPF_NS::random::orientation3D( double& xx, double& yy, double& zz)
//{
//   double x1, x2, sqrt_x1_x2;
//   x1 = uniform_coord(generator);
//   x2 = uniform_coord(generator);
//   
//   while ( x1*x1 + x2*x2 >= 1.0 )
//   {
//      x1 = uniform_coord(generator);
//      x2 = uniform_coord(generator);
//   }
//
//   sqrt_x1_x2 = sqrt(1.0 - x1*x1 - x2*x2);
//   xx = 2.0 * x1 * sqrt_x1_x2;
//   yy = 2.0 * x2 * sqrt_x1_x2;
//   zz = 1.0 - 2.0 * (x1*x1 + x2*x2);
//
//   return EXIT_SUCCESS;
//}

#endif
