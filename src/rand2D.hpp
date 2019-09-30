/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: rand2D.hpp
// Purpose:

#ifndef RAND2D_HPP
#define RAND2D_HPP

#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE
#include <random>    // default_random_engine, uniform_real_distribution

#ifndef PI
#define PI 3.141592653589793238462643382795028814971693993751058209
#endif
#ifndef KB
#define KB 1.380649E-23
#endif

namespace TEM_NS
{
   class random
   {
      public:
         //int rotation( double& omega, double& theta, double& phi);
         //int rotation( double& omega, double& xx, double& yy, double& zz);
         //int orientation3D( double& xx, double& yy, double& zz);
         //int magnitude( double& mm );
         //int displacement1D( double& xx, const double& scale);
         int displacement2D( double& xx, double& yy, const double& scale);
         //int displacement3D( double& xx, double& yy, double& zz, const double& scale);
            //, double& zz)//, const double& T = 1.0)

         std::random_device rd;
         std::mt19937 generator;
         //std::uniform_real_distribution<double> uniform_angle;
         //std::uniform_real_distribution<double> uniform_coord;
         std::normal_distribution<double> normal_distance;

         // CONSTRUCTOR
         random( )
         {
            //std::random_device rd; 
            generator = std::mt19937(rd());
            //uniform_angle
            //   = std::uniform_real_distribution<double>( 0.0, 1.0*PI);
            //uniform_coord
            //   = std::uniform_real_distribution<double>( -1.0, 1.0);
            //uniform_scale
            //   = std::uniform_real_distribution<double>( 0.0, 1.0);
            normal_distance 
               = std::normal_distribution<double>( 0.0, 1.0);
         }
   };
}
#endif
