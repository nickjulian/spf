/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: rand.cpp
// Purpose:

#ifndef RAND_CPP
#define RAND_CPP

#include "rand.hpp"

#include <iostream>
using std::cerr;
using std::endl;

int SPF_NS::random::gaussiandisplacement1D( double& xx, const double& scale)
{
   std::normal_distribution<double> gd( 0.0, 1.0 );
   xx = scale * gd( rr.generator );
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

   xx = uniform_scale( generator ); // uniform real rv in (0,1)

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

//bool SPF_NS::random::bernoulli( const double& pp )
//{
//   std::bernoulli_distribution bb( pp );
//   return bb( generator );
//}

//double SPF_NS::random::exponentialWaitTime( const double& rate, const double& totalRate )
//{
//   double r1;
//   r1 = uniform_scale( generator );
//   return (1.0/totalRate)*log(1.0/r1);
//}

//int SPF_NS::random::rotation( double& omega, 
//                           double& xx, double& yy, double& zz)
//{
//
//   omega = uniform_angle(generator);
//
//   if (! orientation3D( xx, yy, zz) )
//   { 
//      cerr << "Error, SPF_NS::random::orientation3D failed" << endl;
//      return EXIT_FAILURE;
//   }
//
//   return EXIT_SUCCESS;
//}

//int SPF_NS::random::displacement1D( double& xx, const double& scale)
//{
//   xx = scale * normal_distance( generator );
//   return EXIT_SUCCESS;
//}

//int SPF_NS::random::displacement2D( double& xx, double& yy, const double& scale)
//{
//   xx = scale * normal_distance( generator );
//   yy = scale * normal_distance( generator );
//
//   return EXIT_SUCCESS;
//}

//int SPF_NS::random::displacement3D( double& xx, double& yy, double& zz, const double& scale)
//{
//   //double lambda; lambda = 1.0; // damping coefficient 
//
//   //xx = 2 * lambda * KB * T * normal_distance( generator );
//   xx = scale * normal_distance( generator );
//   yy = scale * normal_distance( generator );
//   zz = scale * normal_distance( generator );
//
//   return EXIT_SUCCESS;
//}

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

//bool SPF_NS::random::poisson_trial()
//{
//   if ( uniform_scale(generator) <= subinterval )
//      return true;
//   else
//      return false;
//}
#endif
