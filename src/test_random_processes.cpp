/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include "writeHDF5c.hpp"
#include "rand.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::setw;
using namespace std;
using namespace TEM_NS;

int main( int argc, char* argv[])
{
   size_t Npoints;
   Npoints  = 10;
   
   // increments
   double dt, sqrtdt; dt = 1; sqrtdt = sqrt(dt);
   double dx, dy, dz;
   double rate = 1.0; // Poisson procces rate


   string outFilePrefix;
   outFilePrefix = "gaussian_path";

   //double* odf; odf = new double[Npoints];


   TEM_NS::random rr( rate, dt);

   // Poisson processes
   std::vector<double> Pt;
   //std::vector<double> Px1;
   
   // Wiener processes
   std::vector<double> Wx1; 

   std::vector<double> Wx2;
   std::vector<double> Wy2;

   std::vector<double> Wx3;
   std::vector<double> Wy3;
   std::vector<double> Wz3;

   std::vector<double> Wt;

   // initial coordinate
   Wx1.push_back( 0.0 );

   Wx2.push_back( 0.0 );
   Wy2.push_back( 0.0 );

   Wx3.push_back( 0.0 );
   Wy3.push_back( 0.0 );
   Wz3.push_back( 0.0 );

   Wt.push_back( 0.0 );

   

   for( size_t p=0; p < Npoints; ++p)
   {
      // Poisson process
      if ( rr.poisson_trial() )
      {
         Pt.push_back( p * dt ); 
      }

      // Wiener process
      rr.displacement1D( dx, sqrtdt );
      Wx1.push_back( Wx1.back() + dx );

      rr.displacement2D( dx, dy, sqrtdt );
      Wx2.push_back( Wx2.back() + dx );
      Wy2.push_back( Wy2.back() + dy );

      rr.displacement3D( dx, dy, dz, sqrtdt );
      Wx3.push_back( Wx3.back() + dx );
      Wy3.push_back( Wy3.back() + dy );
      Wz3.push_back( Wz3.back() + dz );

      Wt.push_back( Wt.back() + dt );
   }

   output_path1D_to_hdf5( Wx1, Wt, outFilePrefix );
   output_path2D_to_hdf5( Wx2, Wy2, Wt, outFilePrefix );
   output_path3D_to_hdf5( Wx3, Wy3, Wz3, Wt, outFilePrefix );

   std::cout << "Poisson process event times: " << endl;
   for ( std::vector<double>::iterator Pt_itr = Pt.begin();
         Pt_itr != Pt.end();
         ++Pt_itr)
   {
      cout << *Pt_itr << ", ";
   }
   cout << endl;

   return EXIT_SUCCESS;
}
