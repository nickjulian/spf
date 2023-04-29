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
// File: test_marcus.cpp

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


// integrand specific to this Langevin equation
double integrand(const double& tt, const double& yy, const double& jump)
{
   return jump*yy;
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

   // increments
   double dt; dt = 1;
   //double sqrtdt; sqrtdt = sqrt(dt);

   // rate of the compound poisson process
   double rate; rate = 1.0; 
   double totalRate; totalRate = rate;

   // time and state
   std::vector<double> tt;
   tt.push_back( 0.0 );
   std::vector<double> cc;
   cc.push_back( 1.0 );

   // file will contain times, jump magnitudes, state values
   outFilePrefix = outFilePrefix + "_xdP" + ".h5";

   //TEM_NS::random rr( rate, dt);
   SPF_NS::random rr;

   // a function pointer to pass the integrand to Runge-Kutta
   //double (*zz)(const double& tt, const double& yy, jump);
   //zz = &integrand;

   //int count1, count2; count1 = 0; count2 = 0;
   double wait_time, jumpL, jumpC, jumpRK;

   std::vector<double> unitInterval;
   double rk_dt = 0.01;
   for( size_t t = 0; t < 1.0/rk_dt; ++t)
   {
      unitInterval.push_back( t * rk_dt );
   }

   for( size_t p=0; p < Npoints; ++p)
   {
      // possible implementation for fixed time intervals
      //// Poisson trials, useful if at most one event occurs per interval
      //if ( rr.poisson_trial() )
      //{
      //   tt.push_back( p * dt ); 
      //   //Nmisses = 0;
      //}

      // implementation for variable event time intervals
      wait_time = rr.exponentialWaitTime( rate, totalRate );

      // jump in driving noise
      //rr.displacement1D( jumpL, 1.0 );
      //rr.gaussiandisplacement1D( jumpL, 1.0 );
      rr.unitdisplacement1D( jumpL, 0.5, 0.5 );
      cout << "jumpL: " << jumpL << endl;//debug
      // analytical response in Levy process
      jumpC = jumpDestination( jumpL, cc.back() ) - cc.back();
      cout << "jumpC: " << jumpC << endl;//debug

      // Runge-Kutta response in Levy process
      jumpRK = cc.back(); // initial condition for RK integration
      double rktest; rktest = rk_dt;
      //  loop from 0 to 1 in 
      for( std::vector<double>::const_iterator itr = unitInterval.begin();
            itr != unitInterval.end(); ++itr)
      {
         //cout << "t: " << *itr << endl;//debug
         jumpRK = marcusRK4( integrand, rk_dt, jumpRK, *itr, jumpL);
         rktest = marcusRK4( integrand, rk_dt, rktest, *itr, 1.0);
         //cout << "rktest: " << rktest/rk_dt << endl;
      }
      jumpRK = (jumpRK ) - cc.back();
      cout << "jumpRK: " << jumpRK << endl;//debug
      cout << "rktest: " << rktest << endl << endl;//debug
      tt.push_back( tt.back() + wait_time );
      //cc.push_back( cc.back() + jumpC);
      cc.push_back( cc.back() + jumpRK);
   }

   output_2vectors_to_hdf5( tt, cc, outFilePrefix);
   
   /*****************************************************************/

   //// clear the data, for some reason ... ?
   //tt.clear();
   //cc.clear();

   return EXIT_SUCCESS;
}
