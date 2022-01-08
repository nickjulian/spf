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
// File: writeNetCDF.hpp

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
