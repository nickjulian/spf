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
// File: spf_convergence_to_exponential_1d.cpp

#include <iostream>
#include <fstream>   // ifstream, ofstream
#include <iomanip>
#include <cstdlib>
#include <vector>
#include "writeHDF5c.hpp"
#include "rand1D.hpp"
#include "jumps.hpp"
#include "rungekutta.hpp"
#include "read_cmdline.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::setw;
using namespace std;
using namespace SPF_NS;


// integrand specific to this Langevin equation, b(z)=z , 
//  to use with function pointer.
// In RK4, the integrand is a formula for the time derivative of the 
//  function to extrapolate.
// \Phi_{Y_{k}}(t_{k}^{-}) = \int_{0}^{1} b(z_{0}(\theta))Y_{k} d\theta
double marcusIntegrand( 
      const double& t0, // t0
      const double& y0, // y0
      const double& jump)  // Y_{k}
{
   return jump*y0;
}

//double driftIntegrand(
//               const double& tt, const double& yy, const double& jump)
//{
//   //cout << "tt, yy, jump: " << tt << ", " << yy << ", " << jump << endl;
//   //return (-1.0)*( yy + 0.01*(yy*yy*yy) + cos(yy));
//   return -1.0*yy;
//}

int main( int argc, char* argv[])
{
   // read input parameters
   std::vector<string> args( argv, argv + argc );

   double rate_scale_factor; rate_scale_factor = 1.0;
   double scalar_integrand;
   bool flag_calcstat;
   string output_prefix;
   string input_field_name;
   double dt; 
   int Nt, write_period;
   if (read_cmdline_options(
         args,
         dt,
         Nt,
         scalar_integrand,
         rate_scale_factor,
         write_period,
         flag_calcstat,
         output_prefix,
         input_field_name
            ) != EXIT_SUCCESS )
   {
      cout << "Error: failed to read cmdline arguments" << endl;
      PRINT_USAGE;
      return EXIT_FAILURE;
   }

   // Purpose: to implement a Levy process path_marcus(t) driven by a 
   //          multiplicative Poisson process, (t_k, Y_k) .
   //          L(t) = \sum_{k} Y_{k} U(t-t_{k})    // step functions
   //          path_marcus(t) = path_marcus_{0} 
   //                            + \sum_{k}[ Phi_{k} 
   //                                        - path_marcus(t_{k}^{-}) ]  
   // Loop:
   // 1. determine next jump time time_k and magnitude Y_k, 
   // 2. evaluate the jump response 
   //       Phi_k = z(1) = path_marcus(t_{k}^{-})e^{Y_{k}}
   // 3. increment the Levy process 
   //    path_marcus(t) = path_marcus_{0} 
   //                      + \sum{k}[path_marcus_(t_{k}^{-})e^{Y_{k} 
   //                                  - path_marcus_(t_{k}^{-})]


   double sqrtdt; sqrtdt = sqrt(dt);

   // rate of the compound poisson process
   double w_rate; w_rate = rate_scale_factor;
   double p_rate; p_rate = rate_scale_factor;
   double w_sqrtrate; w_sqrtrate = sqrt(w_rate);

   SPF_NS::random rr;
   std::normal_distribution<double> gd( 0.0, w_sqrtrate );
   std::poisson_distribution<int> pd( dt*p_rate );

   int dP; dP = 0;
   double dB; dB = 0.0;
   double x0; x0 = 1.0;
   // time and state
   std::vector<double> time;
   time.push_back( 0.0 );
   std::vector<double> path_p_marcus;
   std::vector<double> path_p_ito;
   std::vector<double> path_p;   // theoretical solution
   std::vector<double> path_w_ito;
   std::vector<double> path_w_strat;
   std::vector<double> path_w_mil;
   std::vector<double> path_w;
   path_p_marcus.push_back( x0 );
   path_p_ito.push_back( x0 );
   path_p.push_back( x0 );
   path_w_mil.push_back( x0 );
   path_w_strat.push_back( x0 );
   path_w_ito.push_back( x0 );
   path_w.push_back( x0 );
   // processes of driving noises
   double LL; LL = 0.0;//x0;
   double BB; BB = 0.0;//x0; 

   // file will contain times, jump magnitudes, state values
   //output_prefix = output_prefix + "_xdP" + ".h5";

   // a function pointer to pass the integrand to Runge-Kutta
   //double (*zz)(const double& tt, const double& yy, jump);
   //zz = &integrand;

   double jumpC, jumpDestination; // jump magnitudes

   std::vector<double> unitInterval; // domain of marcus mapping integration
   double rk_dt = 0.01;//*dt;
   for( size_t tt = 0; tt < 1.0/rk_dt; ++tt)
   {
      unitInterval.push_back( tt * rk_dt );
   }
   // implementation for fixed time intervals
   for( size_t tt=1; tt < Nt; ++tt) // loop over time steps
   {
      time.push_back( tt*dt);
      // using a Poisson distribution to identify jump count in this dt:
      dP = pd( rr.generator ); 
      
      // update the Poisson process of the driving noise
      //    L(t) = \sum_{k} N(dt,+1)
      LL += dP; 

      // update path_p_ito
      // TODO: double check that this is the ito integral \int X(t)dL(t)
      path_p_ito.push_back( path_p_ito.back() // \sum_{k} x(t_{k})Y_{k} dt
                           + path_p_ito.back() * dP );
                           //+ path_p_ito.back() * dP * dt);//TODO

      // update path_p, being the formula of exact solution 
      //    exp[-0.5 t + L(t)]
      path_p.push_back( x0 * exp( LL ));
      //path_p.push_back( x0 * exp((-0.5* time.back()) + LL ));
      //path_p.append( exp(-0.5* t + path_p_ito.back() ));

      // compute drift component using RK4 over t=t_{n} to t=t_{n} +dt
      //double ddrift;
      // drift = -(x + \epsilon x^3 + cos(x))dt
      //ddrift = marcusRK4( driftIntegrand, // function to integrate
      //                     dt,  // time interval length
      //                     path_p_marcus.back(), // initial position
      //                     time.back(), // time interval beginning
      //                     dP)   // change to position
      //            - path_p_marcus.back();
      //cout << "ddrift: " << ddrift << endl;// debug

      // compute marcus jump response for integrand X(t)Y_{k} dL(t)
      //  by integrating X(t)Y_{k} using Runge Kutta RK4
      jumpDestination = path_p_marcus.back();
      for( std::vector<double>::const_iterator
            itr = unitInterval.begin();
            itr != unitInterval.end(); ++itr)
      {
         jumpDestination = marcusRK4( 
                             marcusIntegrand,//function specific to this SDE
                              rk_dt,   // time increment of unitInterval
                              jumpDestination,  // y_n
                              *itr,    // t_n 
                              dP);  // Y_{k}
      } // jumpDestination now holds y_{n+1} per RK4 notation
      // account for missed drift (makes no difference without drift)
      //jumpDestination = jumpDestination - path_p_marcus.back();
      // If there's drift, add it and jumpDestination and append to 
      // path_p_marcus 
      //  path_p_marcus.append( path_p_marcus.back() 
      //                      + jumpDestination );//+ drift);
      path_p_marcus.push_back( jumpDestination );//+ drift);

      // assign dB = Normal(mean=0,variance=sqrt(w_rate))
      dB = gd( rr.generator );

      // update B(t) = \sum_k \sqrt{dt} dB(t_k)
      BB += sqrtdt * dB;
      cout << "BB: " << BB << endl;//debug

      // update path_w = exp[-0.5*t + B(t)]
      //path_w.append( exp( -0.5 * time.back() + BB) ); // for ito integral
      path_w.push_back( x0 * exp( BB ) );

      // update path_w_ito = \sum_k X(t_k) \sqrt{dt} dB(t_k)
      path_w_ito.push_back( path_w_ito.back() 
                              + path_w_ito.back() * sqrtdt * dB );

      // update path_w_mil
      // Milstein's higher order method (Higham eq. 6.1)
      //
      // W_{next} = W + f(W) dt + g(W) dW + 0.5 g(W) g'(W) ((dW)^2 -dt)
      //  where g(W) = W and f(W) = 0, dW = sqrt(dt)*dB, dB ~ N(0,w_rate)
      // yielding W_t = W_{t-1} + W_{t-1}*dW + 0.5 W_{t-1}((dW)^2 -dt)
      path_w_mil.push_back(
            path_w_mil.back() 
            + ((path_w_mil.back()) * sqrtdt * dB)
            + 0.5*(path_w_mil.back()) * ((dt*dB*dB) - dt)
            );

      // update path_w_strat
      // Ito's lemma yields (Falsone page 204, cites 2 Wong & Zakai papers,
      //  or Kloeden & Platen equation 4.9.7)
      // x = x_0 + \int_0^t 0.5 b'(x) b(x) ds + (I) \int_0^t b(x)dB(s,x)
      //   where b(x) = x, 
      // x = x_0 + (0.5 \sum_s x(s) ds) + \sum_s x(s) sqrt(ds)*dB)
      //   = x_0 + \sum_s x(s)*[ sqrt(dt)*dB + 0.5 dt ]
      //path_w_strat.append(
      //      path_w_strat.back() 
      //      + 0.5 * ((path_w_strat.back())*dt)
      //      + (path_w_strat.back())* sqrtdt * dB
      //      );
      path_w_strat.push_back( 
            path_w_strat.back()
            + (path_w_strat.back()) *( ( 0.5 * dt) + (sqrtdt * dB))
            );
   }

   //output_2vectors_to_hdf5( tt, path_marcus, output_prefix);
   
   /*****************************************************************/
   // open text files to be written to
   ofstream out_text_file_w;
   out_text_file_w.open(output_prefix 
                                 + "_path_w.txt", 
                                 ios::trunc | ios::ate);

   ofstream out_text_file_w_ito;
   out_text_file_w_ito.open(output_prefix 
                                       + "_path_w_ito.txt", 
                                       ios::trunc | ios::ate);

   ofstream out_text_file_w_stratonovich;
   out_text_file_w_stratonovich.open(output_prefix 
                                       + "_path_w_stratonovich.txt", 
                                       ios::trunc | ios::ate);

   ofstream out_text_file_w_milstein;
   out_text_file_w_milstein.open(output_prefix 
                                       + "_path_w_milstein.txt", 
                                       ios::trunc | ios::ate);

   ofstream out_text_file_p_ito;
   out_text_file_p_ito.open(output_prefix 
                                 + "_path_p_ito.txt", 
                                 ios::trunc | ios::ate);

   ofstream out_text_file_p_marcus;
   out_text_file_p_marcus.open(output_prefix 
                                 + "_path_p_marcus.txt", 
                                 ios::trunc | ios::ate);

   ofstream out_text_file_p;
   out_text_file_p.open(output_prefix 
                                 + "_path_p.txt", 
                                 ios::trunc | ios::ate);

   if ( ! out_text_file_w.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix 
               +  "_path_w.txt" 
               << endl;
      out_text_file_w.close();
      return EXIT_FAILURE;
   }

   if ( ! out_text_file_w_ito.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix 
               +  "_path_w_ito.txt" 
               << endl;
      out_text_file_w_ito.close();
      return EXIT_FAILURE;
   }

   if ( ! out_text_file_w_stratonovich.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix 
               +  "_path_w_stratonovich.txt" 
               << endl;
      out_text_file_w_stratonovich.close();
      return EXIT_FAILURE;
   }

   if ( ! out_text_file_w_milstein.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix 
               +  "_path_w_milstein.txt" 
               << endl;
      out_text_file_w_milstein.close();
      return EXIT_FAILURE;
   }

   if ( ! out_text_file_p_ito.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix 
               +  "_path_p_ito.txt" 
               << endl;
      out_text_file_p_ito.close();
      return EXIT_FAILURE;
   }

   if ( ! out_text_file_p.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix 
               +  "_path_p.txt" 
               << endl;
      out_text_file_p.close();
      return EXIT_FAILURE;
   }

   if ( ! out_text_file_p_marcus.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix 
               +  "_path_p_marcus.txt" 
               << endl;
      out_text_file_p_marcus.close();
      return EXIT_FAILURE;
   }
   /*****************************************************************/
   // write the paths to text files
   // loop over time
   if ( time.size() != path_p.size()
         ||
         time.size() != path_p_marcus.size()
         ||
         time.size() != path_p_ito.size()
         ||
         time.size() != path_w.size()
         ||
         time.size() != path_w_ito.size()
         ||
         time.size() != path_w_mil.size()
         ||
         time.size() != path_w_strat.size())
   {
      cout << "Error : time.size() != size of one of paths" << endl;
      return EXIT_FAILURE;
   }

   for( size_t tt=0; tt < time.size(); ++tt)
   {
      out_text_file_w
         << setw(15) << setprecision(8) 
         << time[tt]
         << setw(15) << setprecision(8) 
         << path_w[tt] 
         << endl;

      out_text_file_w_ito
         << setw(15) << setprecision(8) 
         << time[tt]
         << setw(15) << setprecision(8) 
         << path_w_ito[tt] 
         << endl;

      out_text_file_w_stratonovich
         << setw(15) << setprecision(8) 
         << time[tt]
         << setw(15) << setprecision(8) 
         << path_w_strat[tt] 
         << endl;

      out_text_file_w_milstein
         << setw(15) << setprecision(8) 
         << time[tt]
         << setw(15) << setprecision(8) 
         << path_w_mil[tt] 
         << endl;

      out_text_file_p
         << setw(15) << setprecision(8) 
         << time[tt]
         << setw(15) << setprecision(8) 
         << path_p[tt] 
         << endl;

      out_text_file_p_ito
         << setw(15) << setprecision(8) 
         << time[tt]
         << setw(15) << setprecision(8) 
         << path_p_ito[tt] 
         << endl;

      out_text_file_p_marcus
         << setw(15) << setprecision(8) 
         << time[tt]
         << setw(15) << setprecision(8) 
         << path_p_marcus[tt] 
         << endl;
   }

   /*****************************************************************/
   out_text_file_w.close();
   out_text_file_w_ito.close();
   out_text_file_w_stratonovich.close();
   out_text_file_w_milstein.close();
   out_text_file_p.close();
   out_text_file_p_ito.close();
   out_text_file_p_marcus.close();

   /*****************************************************************/

   return EXIT_SUCCESS;
}
