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
// File: rand1D.hpp

#ifndef RAND1D_HPP
#define RAND1D_HPP

#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE
#include <random>    // default_random_engine, uniform_real_distribution

#ifndef PI
#define PI 3.141592653589793238462643382795028814971693993751058209
#endif
//#ifndef KB
//#define KB 1.380649E-23
//#endif

namespace SPF_NS
{
   class random
   {
      public:
         int gaussiandisplacement1D( double& xx, const double& scale);
         int unitdisplacement1D( double& xx, const double& fwd, 
                                    const double& bkwd);
         double exponentialWaitTime( const double& rate, 
                                    const double& totalRate );
         bool bernoulli( const double& pp );
         bool poisson_trial();
         int update_poisson_rate( const double& r, const double& dt )
         {
            //subinterval = 1.0 - exp( -1.0 * r * dt );
            subinterval = r * dt * exp( -1.0 * r * dt );
            return EXIT_SUCCESS;
         }

         std::random_device rd;
         std::mt19937 generator;
         std::uniform_real_distribution<double> uniform_scale;
         std::normal_distribution<double> normal_distance;

         std::poisson_distribution<int> poisson_event_count;

         // CONSTRUCTOR
         random( )
         {
            generator = std::mt19937(rd());
            uniform_scale
               = std::uniform_real_distribution<double>( 0.0, 1.0);
            normal_distance 
               = std::normal_distribution<double>( 0.0, 1.0);
            subinterval = 1.0 * 1.0 * exp( -1.0 * 1.0 * 1.0 );
            poisson_event_count
               = std::poisson_distribution<int>( 1.0 );
         }
         random( const double& r, const double& dt )
         {
            generator = std::mt19937(rd());
            uniform_scale
               = std::uniform_real_distribution<double>( 0.0, 1.0);
            normal_distance 
               = std::normal_distribution<double>( 0.0, 1.0);

            subinterval = r * dt * exp( -1.0 * r * dt );

            poisson_event_count
               = std::poisson_distribution<int>( r );
         }
      private:
         double subinterval;
   };
}
#endif
