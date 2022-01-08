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



int TEM_NS::random::rotation( double& omega, 
                           double& xx, double& yy, double& zz)
{

   omega = uniform_angle(generator);

   if (! orientation3D( xx, yy, zz) )
   { 
      cerr << "Error, TEM_NS::random::orientation3D failed" << endl;
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}

int TEM_NS::random::displacement3D( double& xx, double& yy, double& zz)//, const double& T = 1.0)
{
   //orientation3D( xx, yy, zz);

   double mm;
   double lambda; lambda = 1.0; // damping coefficient 
   //mm = uniform_scale(generator); 

   //xx = 2 * lambda * KB * T * normal_distance( generator );
   xx = normal_distance( generator );
   yy = normal_distance( generator );
   zz = normal_distance( generator );
   //xx = mm * xx;
   //yy = mm * yy;
   //zz = mm * zz;

   return EXIT_SUCCESS;
}


int TEM_NS::random::orientation3D( double& xx, double& yy, double& zz)
{
   double x1, x2, sqrt_x1_x2;
   x1 = uniform_coord(generator);
   x2 = uniform_coord(generator);
   
   while ( x1*x1 + x2*x2 >= 1.0 )
   {
      x1 = uniform_coord(generator);
      x2 = uniform_coord(generator);
   }

   sqrt_x1_x2 = sqrt(1.0 - x1*x1 - x2*x2);
   xx = 2.0 * x1 * sqrt_x1_x2;
   yy = 2.0 * x2 * sqrt_x1_x2;
   zz = 1.0 - 2.0 * (x1*x1 + x2*x2);

   return EXIT_SUCCESS;
}

#endif
