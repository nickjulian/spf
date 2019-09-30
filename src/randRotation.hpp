/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: writeNetCDF.hpp
// Purpose:

#ifndef RANDROTATION_HPP
#define RANDROTATION_HPP

#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE
#include <random>    // default_random_engine, uniform_real_distribution

#ifndef PI
#define PI 3.141592653589793238462643382795028814971693993751058209
#endif

namespace TEM_NS
{
   class random
   {
      public:
         //int rotation( double& omega, double& theta, double& phi);
         int rotation( double& omega, double& xx, double& yy, double& zz);

         //std::random_device rd; // obtains seed for the random number generator
         //std::mt19937 generator(rd()); // mersenne_twister_engine seeded with rd
         //std::uniform_real_distribution<double> uniform_angle(0.0, 2.0*PI);
         //std::uniform_real_distribution<double> uniform_coord(-1.0, 1.0);
         std::random_device rd;
         std::mt19937 generator;
         std::uniform_real_distribution<double> uniform_angle;
         std::uniform_real_distribution<double> uniform_coord;

         // CONSTRUCTOR
         random( )
         {
            //std::random_device rd; 
            generator = std::mt19937(rd());
            uniform_angle
               = std::uniform_real_distribution<double>( 0.0, 1.0*PI);
            uniform_coord
               = std::uniform_real_distribution<double>( -1.0, 1.0);
         }
   };
}
#endif
