/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include "writeHDF5c.hpp"
#include "rand1D.hpp"
#include "jumps.hpp"
#include "rungekutta.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::setw;
using namespace std;
using namespace SPF_NS;


// integrand specific to this Langevin equation, b(z)=z , 
//  to use with function pointer
double integrand(const double& tt, const double& yy, const double& jump)
{
   return jump*yy;
}

double driftIntegrand(
               const double& tt, const double& yy, const double& jump)
{
   //cout << "tt, yy, jump: " << tt << ", " << yy << ", " << jump << endl;
   return (-1.0)*( yy + 0.01*(yy*yy*yy) + cos(yy));
}

int main( int argc, char* argv[])
{
   // read input parameters
   string outFilePrefix;
   std::vector<string> args( argv, argv + argc );
   if ( args.size() != 2 )
   {
      cout << "Usage: " << argv[0] << " <output file prefix>" << endl;
      return EXIT_FAILURE;
   }
   else if ( args.size() == 2 )
   {
      outFilePrefix = string( args[1] );
   }

   // Purpose: to implement a Levy process c(t) driven by a multiplicative 
   //          Poisson process, (t_k, Y_k) .
   //          L(t) = \sum_{k} Y_{k} U(t-t_{k})
   //          c(t) = c_{0} + \sum_{k}[ Phi_{k} - c(t_{k}^{-}) ]
   // Loop:
   // 1. determine next jump time t_k and magnitude Y_k, 
   // 2. evaluate the jump response Phi_k = z(1) = c(t_{k}^{-})e^{Y_{k}}
   // 3. increment the Levy process 
   //    c(t) = c_{0} + \sum{k}[c(t_{k}^{-})e^{Y_{k} - c(t_{k}^{-})]

   int Npoints;
   Npoints  = 100;
   double runTime = 2.0;

   // increments
   double dt; dt = runTime/Npoints;//0.01;
   //double sqrtdt; sqrtdt = sqrt(dt);

   // rate of the compound poisson process
   double rate; rate = 2.5*dt;// 0.01; 
   double totalRate; totalRate = rate;

   // time and state
   std::vector<double> tt;
   tt.push_back( 0.0 );
   std::vector<double> cc;
   cc.push_back( 1.0 );

   // file will contain times, jump magnitudes, state values
   outFilePrefix = outFilePrefix + "_xdP" + ".h5";

   SPF_NS::random rr( rate, dt);
   //SPF_NS::random rr;

   // a function pointer to pass the integrand to Runge-Kutta
   //double (*zz)(const double& tt, const double& yy, jump);
   //zz = &integrand;

   //int count1, count2; count1 = 0; count2 = 0;
   double wait_time, jumpL, jumpC, jumpRK;
   int jumpCount;
   double jumpTotal; 
   std::vector<double> unitInterval;
   double rk_dt = 0.01;
   for( size_t t = 0; t < 1.0/rk_dt; ++t)
   {
      unitInterval.push_back( t * rk_dt );
   }

   for( size_t p=0; p < Npoints; ++p)
   {
      //// implementation for fixed time intervals
      // using a Poisson distribution to identify jump count in this dt
      jumpCount = rr.poisson_event_count( rr.generator );
      cout << "jumpCount: " << jumpCount << endl;//debug

      // compute drift component
      double ddrift;
      // drift = -(x + \epsilon x^3 + cos(x))dt
      ddrift = marcusRK4( driftIntegrand, dt, cc.back(), tt.back(), 1.0) 
                           - cc.back();
      cout << "ddrift: " << ddrift << endl;// debug

      jumpTotal = 0.0;
      // compute jump component
      for( size_t m=0; m < jumpCount; ++m)
      {
         // jump in driving noise
         //rr.gaussiandisplacement1D( jumpL, 1.0 );
         rr.unitdisplacement1D( jumpL, 0.5, 0.5 );
         //cout << "jumpL: " << jumpL << endl;//debug
         // analytical response in Levy process
         //jumpC = jumpDestination( jumpL, cc.back() ) - cc.back();
         //cout << "jumpC: " << jumpC << endl;//debug

         // Runge-Kutta response in Levy process
         jumpRK = cc.back(); // initial condition for RK integration
         //double rktest; rktest = rk_dt;
         // integrate the marcus mapping 
         for( std::vector<double>::const_iterator 
                  itr = unitInterval.begin(); 
               itr != unitInterval.end(); ++itr)
         {
            //cout << "t: " << *itr << endl;//debug
            // iterate over unit interval increments
            jumpRK = marcusRK4( integrand, 
                                 rk_dt, // time increment of unit interval
                                 jumpRK,  // y_n for rk4 
                                 *itr,    // t_n for rk4
                                 jumpL);  // Y_k
            //rktest = marcusRK4( integrand, rk_dt, rktest, *itr, 1.0);
            //cout << "rktest: " << rktest/rk_dt << endl;
         }
         jumpRK = (jumpRK ) - cc.back();
         jumpTotal = jumpTotal + jumpRK;
         //tt.push_back( p*dt);
         //cc.push_back( cc.back() + jumpRK);
         cout << "jumpRK: " << jumpRK << endl;//debug
      }
      tt.push_back( p*dt);
      cc.push_back( cc.back() + jumpTotal + ddrift);
      // add the drift part
      //tt.push_back( tt.back() + dt);
      //cc.push_back( cc.back() + ddrift );
   }

   output_2vectors_to_hdf5( tt, cc, outFilePrefix);
   
   /*****************************************************************/

   //// clear the data, for some reason ... ?
   //tt.clear();
   //cc.clear();

   return EXIT_SUCCESS;
}
