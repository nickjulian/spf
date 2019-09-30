/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: writeNetCDF.cpp
// Purpose:

#ifndef RANDROTATION_CPP
#define RANDROTATION_CPP

#include "randRotation.hpp"



int TEM_NS::random::rotation( double& omega, 
                           double& xx, double& yy, double& zz)
{
   //std::random_device rd; // obtaines seed for the random number generator
   //std::mt19937 generator(rd()); // mersenne_twister_engine seeded with rd
   //std::uniform_real_distribution<double> uniform_angle(0.0, 2.0*PI);
   //std::uniform_real_distribution<double> uniform_coord(-1.0, 1.0);

   double x1, x2, sqrt_x1_x2;
   omega = uniform_angle(generator);
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
