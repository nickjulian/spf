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
// File: test_random_processes.cpp

#include <iostream>  // cout, cin, cerr, endl
#include <iomanip>   // setw, setprecision
#include <fstream>   // ifstream, ofstream
#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <string>
#include <vector>
#include "writeHDF5c.hpp"
#include "rand.hpp"
#include "read_cmdline.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::setw;
using namespace std;
using namespace SPF_NS;

int main( int argc, char* argv[])
{
   int Nt;
   double rate_scale_factor;
   double dt;
   int write_period;
   bool flag_calcstat;
   string output_prefix;
   string input_field_name;
   
   // increments
   double dt, sqrtdt; 
   int Nt, write_period
   //double dx, dy, dz;
   double rate = 1.0; // Poisson procces rate


   string outFilePrefix;
   outFilePrefix = "gaussian_path";
   std::vector<string> args( argv, argv + argc );

   if (read_cmdline_options(
         args,
         dt,
         Nt,
         rate_scale_factor,
         write_period,
         flag_calcstat,
         output_prefix,
         input_field_name
            ) != EXIT_SUCCESS )
   {
      cout << "Error: failed to read cmdline arguments" << endl;
      return EXIT_FAILURE;
   }
   sqrtdt = sqrt(dt);

   std::vector<double> Tt; // time 

   double rate; rate = rate_scale_factor;
   double drift; drift = rate;

   // random generator
   SPF_NS::random rr;

   // Poisson processes
   int dP;
   std::poisson_distribution<int> pd( rate );
   //std::vector<double> Pt;
   std::list<std::vector<double>> Pp;
   //std::list<std::vector<double>> Xp; // dX(t) is the SDE driven by dP(t)
   //std::list<std::vector<double>> Xp_exact;// path according to analysis
   
   // Wiener processes
   double dW;
   std::normal_distribution<double> gd( 0.0, 1.0);
   double dW_withdrift;
   double dW_withdrift2;
   std::normal_distribution<double> 
      gd_withdrift2( dt * drift, sqrtdt * rate );

   std::list<std::vector<double>> Ww; // W(t) is the stochastic process
   std::list<std::vector<double>> Ww_withdrift;
   //std::list<std::vector<double>> Xw; // dX(t) is the SDE driven by dW(t)
   //std::list<std::vector<double>> Xw_exact;// path according to analysis

   //std::vector<double>::iterator Tt_itr;
   std::list<std::vector<double>>::iterator Ww_itr;
   std::list<std::vector<double>>::iterator Ww_withdrift_itr;
   std::list<std::vector<double>>::iterator Pp_itr;
   //std::list<std::vector<double>>::iterator Xx_itr;
   //std::list<std::vector<double>>::iterator Xx_exact_itr;

   Tt.push_back( 0.0 );
   Ww.push_back( new std::vector<double> );
   Ww_withdrift.push_back( new std::vector<double> );
   Pp.push_back( new std::vector<double> );
   //Xx.push_back( new std::vector<double> );
   //Xx_exact.push_back( new std::vector<double> );

   Ww_itr = Ww.begin();
   Ww_withdrift_itr = Ww_withdrift.begin();
   Pp_itr = Pp.begin();
   //Xx_itr = Xx.begin();
   //Xx_exact_itr = Xx_exact.begin();

   // initial coordinate
   Ww_itr->push_back( 0.0 );
   Ww_withdrift_itr->push_back( 0.0 );
   Pp_itr->push_back( 0.0 );
   //Xx_itr->push_back( 0.0 );
   //Xx_exact_itr->push_back( 0.0 );

   // text file output
   ofstream out_text_file_w;
   out_text_file_w.open(outFilePrefix + "_path_w.txt", ios:app | ios::ate);

   ofstream out_text_file_w_withdrift;
   out_text_file_w_withdrift.open(
         outFilePrefix + "_path_w_withdrift.txt", ios:app | ios::ate);

   ofstream out_text_file_p;
   out_text_file_p.open(outFilePrefix + "_path_p.txt", ios:app | ios::ate);

   if ( ! out_text_file_p.good() )
   {
      cout << "warning: could not open output file: "
            << outFilePrefix + "_path_p.txt" << endl;
      out_text_file_p.close();
      return EXIT_FAILURE;
   }
   if ( ! out_text_file_w.good() )
   {
      cout << "warning: could not open output file: "
            << outFilePrefix + "_path_w.txt" << endl;
      out_text_file_w.close();
      return EXIT_FAILURE;
   }

   for( size_t tt=0; tt < Nt; ++tt)
   {
      Tt.push_back( Tt.back() + dt );
      // Poisson process
      //if ( rr.poisson_trial() )
      //{
      //   Pt.push_back( p * dt ); 
      //}
      dP = pd( rr.generator );

      Pp_itr->push_back( dP );

      // Wiener process
      // dW ~ \sqrt(\Delta t) * N(0, 1.0)
      // W_t = W_{t-1} + f(x) dt + g(x) dW + 0.5 g(x) g'(x) ((dW)^2 -dt)
      // g(x) = rate_scale_factor * ?something, which here is 1?
      // g'(x) = d/dx g(x) = 0
      // f(x) = drift
      dW = gd( rr.generator );
      Ww_itr->push_back( 
               Ww_itr->back() +  rate * sqrtdt * dW // + 0, since g'(x) = 0
            );

      // dW ~ \sqrt(\Delta t) * N( 0, 1.0) 
      dW_withdrift = gd( rr.generator );
      Ww_withdrift_itr->push_back( 
               Ww_withdrift_itr->back() 
                  + (drift * dt) + (rate * sqrtdt * dW_withdrift )
            );

      // dW ~ \sqrt(\Delta t) * N( drift, rate) // may be wrong w/drift
      dW_withdrift2 = gd_withdrift2( rr.generator );
      Ww_withdrift2_itr->push_back( 
               Ww_withdrift_itr2->back() 
                  + (dW_withdrift2 ) // drift and rate already in dW
            );

      // TODO: evaluate the analytical solution
      // The SDE
      // dX(t) = a X(t) dt + \mu X(t) dW(t)
      // has exact solution
      // X(t) = X(0) exp( (a - 0.5 \mu^2)t + \mu W(t))
      
      // append the positions and times to text files
      out_text_file_w_withdrift  
         << setw(12) << setprecision(8)
         *(Tt.back())
         << setw(12) << setprecision(8)
         *(Pp_itr->back())
         << setw(12) << setprecision(8)
         *(Ww_itr->back())
         << setw(12) << setprecision(8)
         *(Ww_itr->back())
         << setw(12) << setprecision(8)
         *(Ww_withdrift_itr->back())
         << setw(12) << setprecision(8)
         *(Ww_withdrift2_itr->back())
         << endl;
   }

   //output_path1D_to_hdf5( Wx1, Wt, outFilePrefix );
   //output_path2D_to_hdf5( Wx2, Wy2, Wt, outFilePrefix );
   //output_path3D_to_hdf5( Wx3, Wy3, Wz3, Wt, outFilePrefix );

   //std::cout << "Poisson process event times: " << endl;
   //for ( std::vector<double>::iterator Pt_itr = Pt.begin();
   //      Pt_itr != Pt.end();
   //      ++Pt_itr)
   //{
   //   cout << *Pt_itr << ", ";
   //}
   //cout << endl;

   out_text_file_w.close();
   out_text_file_p.close();

   for (size_t ii= 0; ii < Ww.size(); ++ii) 
   {
      delete *(Ww.back());
      Ww.pop_back();
   }
   for (size_t ii= 0; ii < Ww_withdrift.size(); ++ii) 
   {
      delete *(Ww_withdrift.back());
      Ww_withdrift.pop_back();
   }
   for (size_t ii= 0; ii < Pp.size(); ++ii) 
   {
      delete *(Pp.back());
      Pp.pop_back();
   }
   //for (size_t ii= 0; ii < Xx.size(); ++ii) 
   //{
   //   delete *(Xx.back());
   //   Xx.pop_back();
   //}
   //for (size_t ii= 0; ii < Xx_exact.size(); ++ii) 
   //{
   //   delete *(Xx_exact.back());
   //   Xx_exact.pop_back();
   //}
   return EXIT_SUCCESS;
}
