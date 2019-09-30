/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: rand2D.cpp
// Purpose:

#ifndef RAND2D_CPP
#define RAND2D_CPP

#include "rand2D.hpp"


//int TEM_NS::random::rotation( double& omega, 
//                           double& xx, double& yy, double& zz)
//{
//
//   omega = uniform_angle(generator);
//
//   if (! orientation3D( xx, yy, zz) )
//   { 
//      cerr << "Error, TEM_NS::random::orientation3D failed" << endl;
//      return EXIT_FAILURE;
//   }
//
//   return EXIT_SUCCESS;
//}

int TEM_NS::random::displacement2D( double& xx, double& yy, const double& scale);//, double& zz)//, const double& T = 1.0)
{
   //orientation3D( xx, yy, zz);

   //double mm;
   //double lambda; lambda = 1.0; // damping coefficient 
   //mm = uniform_scale(generator); 

   //xx = 2 * lambda * KB * T * normal_distance( generator );
   xx = scale * normal_distance( generator );
   yy = scale * normal_distance( generator );
   //zz = normal_distance( generator );
   //xx = mm * xx;
   //yy = mm * yy;
   //zz = mm * zz;

   return EXIT_SUCCESS;
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
