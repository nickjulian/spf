/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: rand1D.cpp
// Purpose:

#ifndef RAND1D_CPP
#define RAND1D_CPP

#include "rand1D.hpp"

#include <iostream>
using std::cout;
using std::endl;

double TEM_NS::random::exponentialWaitTime( const double& rate, const double& totalRate )
{
   double r1;
   r1 = uniform_scale( generator );
   return (1.0/totalRate)*log(1.0/r1);
}

bool TEM_NS::random::poisson_trial()
{
   if ( uniform_scale(generator) <= subinterval )
      return true;
   else
      return false;
}

int TEM_NS::random::gaussiandisplacement1D( double& xx, const double& scale)
{
   xx = scale * (2.0*uniform_scale( generator ) - 1.0);
   return EXIT_SUCCESS;
}

int TEM_NS::random::unitdisplacement1D( double& xx, const double& fwd, 
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

bool TEM_NS::random::bernoulli( const double& pp )
{
   std::bernoulli_distribution bb( pp );
   return bb( generator );
}

//int TEM_NS::random::orientation3D( double& xx, double& yy, double& zz)
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
