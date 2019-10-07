/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: rand.hpp
// Purpose:

#ifndef RAND_HPP
#define RAND_HPP

#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE
#include <random>    // default_random_engine, uniform_real_distribution

#ifndef PI
#define PI 3.141592653589793238462643382795028814971693993751058209
#endif
#ifndef KB
#define KB 1.380649E-23
#endif

namespace SPF_NS
{
   class random
   {
      public:
         //int rotation( double& omega, double& theta, double& phi);
         //int rotation( double& omega, double& xx, double& yy, double& zz);
         //int orientation3D( double& xx, double& yy, double& zz);
         ////int magnitude( double& mm );
         //int displacement1D( double& xx, const double& scale);
         //int displacement2D( double& xx, double& yy, const double& scale);
         //int displacement3D( double& xx, double& yy, double& zz, const double& scale);
            //, double& zz)//, const double& T = 1.0)
         
         int gaussiandisplacement1D( double& xx, const double& scale);
         int unitdisplacement1D( double& xx, const double& fwd, 
                                    const double& bkwd);
         //double exponentialWaitTime( const double& rate, 
         //                           const double& totalRate );
         //bool bernoulli( const double& pp );

         // poisson process 
         //bool poisson_trial();
         
         //double get_subinterval( )
         //{
         //   return subinterval;
         //}

         //int update_poisson_rate( const double& r, const double& dt )
         //{
         //   // free the previous poisson_event_count and create a new one
         //   poisson_event_count
         //      = std::poisson_distribution<int>( r );
         //   ////subinterval = 1.0 - exp( -1.0 * r * dt );
         //   //subinterval = r * dt * exp( -1.0 * r * dt );
         //   return EXIT_SUCCESS;
         //}

         std::random_device rd;
         std::mt19937 generator;
         ////std::uniform_real_distribution<double> 
         ////   uniform_positive_unit_distribution;
         ////std::uniform_real_distribution<double> uniform_angle;
         ////std::uniform_real_distribution<double> uniform_coord;
         std::uniform_real_distribution<double> uniform_scale;
         ////std::normal_distribution<double> normal_distance;
         //std::normal_distribution<double> gaussian_sample;
         //std::poisson_distribution<int> poisson_event_count;
         // CONSTRUCTOR
         random( const double& r, const double& dt )
         {
            //std::random_device rd; 
            generator = std::mt19937(rd());
            //poisson_event_count
            //   = std::poisson_distribution<int>( r );
            //gaussian_sample
            //   = std::normal_distribution<double>( 0.0, r );
            //uniform_angle
            //   = std::uniform_real_distribution<double>( 0.0, 1.0*PI);
            //uniform_coord
            //   = std::uniform_real_distribution<double>( -1.0, 1.0);
            //normal_distance 
            //   = std::normal_distribution<double>( 0.0, 1.0);
            //uniform_positive_unit_distribution
            //   = std::uniform_real_distribution<double>( 0.0, 1.0);
            uniform_scale
               = std::uniform_real_distribution<double>( 0.0, 1.0);

            ////subinterval = 1.0 - exp( -1.0 * r * dt );
            //subinterval = r * dt * exp( -1.0 * r * dt );
         }

      private:
         //double subinterval;
   };
}
#endif
