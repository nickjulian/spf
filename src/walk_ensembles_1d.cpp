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
// File: walk_ensembles_1d.cpp

#include <iostream>  // cout, cin, cerr, endl
#include <iomanip>   // setw, setprecision
#include <fstream>   // ifstream, ofstream
#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <string>
#include <vector>
#include <list>
#include "writeHDF5c.hpp"
#include "rand.hpp"
#include "read_cmdline.hpp"

//#define PRINT_USAGE cout << "Usage: " << argv[0] << " <OPTIONS>" << endl << "OPTIONS : " << endl << "   -i <input field hdf5 file>" << endl << "   -o <output file prefix>" << endl << "   -Nt <number of time steps>" << endl << "   [-dt <time increment>]" << endl << "   [-r <rate scale factor>]" << endl << "   [-wp <steps between file writes>]" << endl << "   [-stat]"<< endl;

using std::cout;
using std::endl;
using std::cerr;
using std::setw;
using std::setprecision;
using namespace std;
using namespace SPF_NS;

int main( int argc, char* argv[])
{
   double rate_scale_factor; rate_scale_factor = 1.0;
   double scalar_integrand;
   bool flag_calcstat;
   string output_prefix;
   string input_field_name;
   
   // increments
   double dt, sqrtdt; 
   //size_t Nensemble; Nensemble = 100000;
   size_t Nensemble; Nensemble = 1;
   int Nt, write_period;
   //double dx, dy, dz;

   std::vector<string> args( argv, argv + argc );

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

   std::vector<double> Tt; // time 
   double rate; rate = rate_scale_factor;
   double sqrtrate; 
   double drift; drift = rate * scalar_integrand;

   sqrtdt = sqrt(dt);
   sqrtrate = sqrt(rate);

   // random generator
   SPF_NS::random rr;

   // Poisson processes
   int dP;
   std::poisson_distribution<int> pd( rate );
   //std::vector<double> Pt;
   std::list<std::vector<double>*> Pp;
   //std::list<std::vector<double>*> Xp; // dX(t) is the SDE driven by dP(t)
   //std::list<std::vector<double>*> Xp_exact;// path according to analysis
   
   // Wiener processes
   double dW;
   std::normal_distribution<double> gd( 0.0, sqrtrate); // (mean, stddev)
   double dW_withdrift;
   double dW_withdrift2;
   std::normal_distribution<double> 
      gd_withdrift2( dt * drift, sqrtdt * sqrtrate );

   std::list<std::vector<double>*> Ww; // W(t) is the stochastic process
   std::list<std::vector<double>*> Ww_withdrift;
   std::list<std::vector<double>*> Ww_withdrift2;
   //std::list<std::vector<double>*> Xw; // dX(t) is the SDE driven by dW(t)
   //std::list<std::vector<double>*> Xw_exact;// path according to analysis

   //std::vector<double>::iterator Tt_itr;
   std::vector<double>* Ww_ptr;
   std::vector<double>* Ww_withdrift_ptr;
   std::vector<double>* Ww_withdrift2_ptr;
   std::vector<double>* Pp_ptr;
   //std::vector<double> Xx_ptr;
   //std::vector<double> Xx_exact_ptr;

   // text file output
   ofstream out_text_file_w;
   out_text_file_w.open(output_prefix + "_path_w.txt", ios::app | ios::ate);

   ofstream out_text_file_w_withdrift;
   out_text_file_w_withdrift.open(
         output_prefix + "_path_w_withdrift.txt", ios::app | ios::ate);

   ofstream out_text_file_w_withdrift2;
   out_text_file_w_withdrift2.open(
         output_prefix + "_path_w_withdrift2.txt", ios::app | ios::ate);

   ofstream out_text_file_p;
   out_text_file_p.open(output_prefix + "_path_p.txt", ios::app | ios::ate);

   ofstream out_text_file_p_stats;
   out_text_file_p_stats.open(
         output_prefix + "_path_p_stats.txt", ios::app | ios::ate);

   ofstream out_text_file_w_stats;
   out_text_file_w_stats.open(
         output_prefix + "_path_w_stats.txt", ios::app | ios::ate);

   ofstream out_text_file_w_withdrift_stats;
   out_text_file_w_withdrift_stats.open(
         output_prefix + "_path_w_withdrift_stats.txt", 
         ios::app | ios::ate);

   ofstream out_text_file_w_withdrift2_stats;
   out_text_file_w_withdrift2_stats.open(
         output_prefix + "_path_w_withdrift2_stats.txt", 
         ios::app | ios::ate);

   if ( ! out_text_file_p.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix + "_path_p.txt" << endl;
      out_text_file_p.close();
      return EXIT_FAILURE;
   }
   if ( ! out_text_file_w.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix + "_path_w.txt" << endl;
      out_text_file_w.close();
      return EXIT_FAILURE;
   }
   if ( ! out_text_file_w_withdrift.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix + "_path_w_withdrift.txt" << endl;
      out_text_file_w_withdrift.close();
      return EXIT_FAILURE;
   }
   if ( ! out_text_file_w_withdrift2.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix + "_path_w_withdrift2.txt" << endl;
      out_text_file_w_withdrift2.close();
      return EXIT_FAILURE;
   }

   cout << "(scalar_integrand, rate) : (" << scalar_integrand << ", " << rate << ")" << endl; // debug

   Tt.push_back( 0.0 );
   for( int tt=0; tt < Nt; ++tt) 
   {
      if ( tt % write_period == 0 ) 
      {
         Tt.push_back( Tt[Tt.size() -1] + dt );
      }
      else Tt[Tt.size()-1] = Tt[Tt.size()-1] + dt;
   }
   cout << "Tt.size() : " << Tt.size() << endl; // debug

   // loop over ensemble members
   for ( size_t en_count=0; en_count < Nensemble; ++en_count)
   {
      Ww.push_back( new std::vector<double> );
      Ww_withdrift.push_back( new std::vector<double> );
      Ww_withdrift2.push_back( new std::vector<double> );
      Pp.push_back( new std::vector<double> );
      //Xx.push_back( new std::vector<double> );
      //Xx_exact.push_back( new std::vector<double> );

      Ww_ptr = Ww.back();
      Ww_withdrift_ptr = Ww_withdrift.back();
      Ww_withdrift2_ptr = Ww_withdrift2.back();
      Pp_ptr = Pp.back();
      //Xx_ptr = Xx.begin();
      //Xx_exact_ptr = Xx_exact.begin();

      // initial coordinate
      Ww_ptr->push_back( 0.0 );
      Ww_withdrift_ptr->push_back( 0.0 );
      Ww_withdrift2_ptr->push_back( 0.0 );
      Pp_ptr->push_back( 0.0 );
      //Xx_ptr->push_back( 0.0 );
      //Xx_exact_ptr->push_back( 0.0 );

      // loop over time
      for( int tt=0; tt < Nt; ++tt)
      {
         // TODO: evaluate the analytical solution
         // The SDE
         // dX(t) = a X(t) dt + \mu X(t) dW(t)
         // has exact solution
         // X(t) = X(0) exp( (a - 0.5 \mu^2)t + \mu W(t))
         
         // Poisson process
         //if ( rr.poisson_trial() )
         //{
         //   Pt.push_back( p * dt ); 
         //}
         dP = pd( rr.generator );

         if ( tt % write_period == 0 )
         {
            Pp_ptr->push_back( (*Pp_ptr)[ Pp_ptr->size() -1]
                                 + (scalar_integrand * dt * dP) );
         }
         else
         {
            (*Pp_ptr)[Pp_ptr->size() -1] = (*Pp_ptr)[Pp_ptr->size() -1] 
                                             + (scalar_integrand * dt * dP);
         }

         // Wiener process
         // dW ~ \sqrt(\Delta t) * N(0, 1.0)
         // W_t = W_{t-1} + f(x) dt + g(x) dW + 0.5 g(x) g'(x) ((dW)^2 -dt)
         // g(x) = rate_scale_factor * ?something, which here is 1?
         // g'(x) = d/dx g(x) = 0
         // f(x) = drift
         dW = gd( rr.generator );
         if ( tt % write_period == 0 )
         {
            Ww_ptr->push_back( 
                  (*Ww_ptr)[Ww_ptr->size() -1]
                     + (scalar_integrand * sqrtdt 
                         * dW) // + 0, since g'(x) =0
                         //* sqrtrate * dW) // + 0, since g'(x) =0
               );
         }
         else
         {
            (*Ww_ptr)[Ww_ptr->size() -1]
                  = (*Ww_ptr)[Ww_ptr->size() -1]
                              + (scalar_integrand * sqrtdt 
                                  * dW); // + 0, since g'(x) =0
         }

         // dW ~ \sqrt(\Delta t) * N( 0, 1.0) 
         dW_withdrift = gd( rr.generator );
         if ( tt % write_period == 0 )
         {
            Ww_withdrift_ptr->push_back( 
                  (*Ww_withdrift_ptr)[Ww_withdrift_ptr->size() -1]
                     + (drift * dt) 
                     + (scalar_integrand * sqrtdt * dW_withdrift)
                     //+ (scalar_integrand * sqrtdt * sqrtrate * dW_withdrift)
               );
         }
         else
         {
            (*Ww_withdrift_ptr)[Ww_withdrift_ptr->size() -1]
               = (*Ww_withdrift_ptr)[Ww_withdrift_ptr->size() -1]
                     + (drift * dt) 
                     + (scalar_integrand * sqrtdt * dW_withdrift);
         }

         // dW ~ \sqrt(\Delta t) * N( drift, rate) // may be wrong w/drift
         dW_withdrift2 = gd_withdrift2( rr.generator );
         if ( tt % write_period == 0 )
         {
            Ww_withdrift2_ptr->push_back( 
                  (*Ww_withdrift2_ptr)[Ww_withdrift2_ptr->size() -1]
                     + (scalar_integrand * dW_withdrift2 ) // drift and rate already in dW
               );
         }
         else
         {
            (*Ww_withdrift2_ptr)[Ww_withdrift2_ptr->size() -1]
               = (*Ww_withdrift2_ptr)[Ww_withdrift2_ptr->size() -1]
                     + (scalar_integrand * dW_withdrift2 ); // drift and rate already in dW
         }

         Pp_ptr = Pp.back();
         Ww_ptr = Ww.back();
         Ww_withdrift_ptr = Ww_withdrift.back();
         Ww_withdrift2_ptr = Ww_withdrift2.back();
      } // loop over time
      //cout << "path # " << en_count << ", Pp_ptr->size() : " << Pp_ptr->size() << endl; // debug
   } // loop over ensemble members

   // write the results to files and calculate moments
   double pp_sum, pp_sqrsum;
   double ww_sum, ww_sqrsum; 
   double ww_withdrift_sum, ww_withdrift_sqrsum;
   double ww_withdrift2_sum, ww_withdrift2_sqrsum;

   std::vector<double>::const_iterator Tt_itr;
   Tt_itr = Tt.begin();
   for( int tt=0; tt < Tt.size(); ++tt)
   //for ( int tt=0; tt < Nt; ++tt )
   {
      pp_sum = 0.0; pp_sqrsum = 0.0;
      ww_sum = 0.0; ww_sqrsum = 0.0;
      ww_withdrift_sum = 0.0; ww_withdrift_sqrsum = 0.0;
      ww_withdrift2_sum = 0.0; ww_withdrift2_sqrsum = 0.0;

      // append the positions and times to text files
      //if ( tt % write_period == 0 )
      //{
         out_text_file_p
            << setw(12) << setprecision(8)
            << *Tt_itr;
      //}
      // increment through ensemble members
      for (std::list<std::vector<double>*>::iterator 
            itr = Pp.begin(); itr != Pp.end(); ++itr)
      {
         Pp_ptr = *itr;
         //if ( tt % write_period == 0 )
         //{
            out_text_file_p << setw(12) << setprecision(8)
               << (*Pp_ptr)[tt];
         //}
         pp_sum += (*Pp_ptr)[tt];
         pp_sqrsum += ((*Pp_ptr)[tt]) * ((*Pp_ptr)[tt]);
      }

      //if ( tt % write_period == 0 )
      //{
         out_text_file_w
            << setw(12) << setprecision(5)
            << *Tt_itr;
      //}
      for (std::list<std::vector<double>*>::iterator 
            itr = Ww.begin(); itr != Ww.end(); ++itr)
      {
         Ww_ptr = *itr;
         //if ( tt % write_period == 0 )
         //{
            out_text_file_w << setw(12) << setprecision(5)
               << (*Ww_ptr)[tt];
         //}
         ww_sum += (*Ww_ptr)[tt];
         ww_sqrsum += (*Ww_ptr)[tt] * (*Ww_ptr)[tt];
      }

      //if ( tt % write_period == 0 )
      //{
         out_text_file_w_withdrift
            << setw(12) << setprecision(5)
            << *Tt_itr;
      //}
      for (std::list<std::vector<double>*>::iterator 
            itr = Ww_withdrift.begin(); itr != Ww_withdrift.end(); ++itr)
      {
         Ww_ptr = *itr;
         //if ( tt % write_period == 0 )
         //{
            out_text_file_w_withdrift << setw(12) << setprecision(5)
               << (*Ww_ptr)[tt];
         //}
         ww_withdrift_sum += (*Ww_ptr)[tt];
         ww_withdrift_sqrsum += (*Ww_ptr)[tt] * (*Ww_ptr)[tt];
      }

      //if ( tt % write_period == 0 )
      //{
         out_text_file_w_withdrift2
            << setw(12) << setprecision(5)
            << *Tt_itr;
      //}
      for (std::list<std::vector<double>*>::iterator 
            itr = Ww_withdrift2.begin(); itr != Ww_withdrift2.end(); ++itr)
      {
         Ww_ptr = *itr;
         //if ( tt % write_period == 0 )
         //{
            out_text_file_w_withdrift2 << setw(12) << setprecision(5)
               << (*Ww_ptr)[tt];
         //}
         ww_withdrift2_sum += (*Ww_ptr)[tt];
         ww_withdrift2_sqrsum += (*Ww_ptr)[tt] * (*Ww_ptr)[tt];
      }

      //if ( tt % write_period == 0 )
      //{
         out_text_file_p << endl;
         out_text_file_w << endl;
         out_text_file_w_withdrift << endl;
         out_text_file_w_withdrift2 << endl;

         out_text_file_p_stats
            << setw(12) << setprecision(5)
            << *Tt_itr;
         out_text_file_p_stats 
            << setw(12) << setprecision(5)
            << pp_sum/Nensemble 
            << setw(12) << setprecision(5)
            << (pp_sqrsum/Nensemble) 
                  - ((pp_sum/Nensemble) * (pp_sum/Nensemble))
            << endl;

         out_text_file_w_stats
            << setw(12) << setprecision(5)
            << *Tt_itr;
         out_text_file_w_stats 
            << setw(12) << setprecision(5)
            << ww_sum/Nensemble 
            << setw(12) << setprecision(5)
            << (ww_sqrsum/Nensemble) 
                  - ((ww_sum/Nensemble) * (ww_sum/Nensemble))
            << endl;

         out_text_file_w_withdrift_stats
            << setw(12) << setprecision(5)
            << *Tt_itr;
         out_text_file_w_withdrift_stats 
            << setw(12) << setprecision(5)
            << ww_withdrift_sum/Nensemble 
            << setw(12) << setprecision(5)
            << (ww_withdrift_sqrsum/Nensemble)
               - ((ww_withdrift_sum/Nensemble) 
                     * (ww_withdrift_sum/Nensemble))
            << endl;

         out_text_file_w_withdrift2_stats
            << setw(12) << setprecision(5)
            << *Tt_itr;
         out_text_file_w_withdrift2_stats 
            << setw(12) << setprecision(5)
            << ww_withdrift2_sum/Nensemble 
            << setw(12) << setprecision(5)
            << (ww_withdrift2_sqrsum/Nensemble)
               - ((ww_withdrift2_sum/Nensemble) 
                     * (ww_withdrift2_sum/Nensemble))
            << endl;
      //}

      //if ( tt % write_period == 0 ) 
         ++Tt_itr;
   }

   //std::cout << "Poisson process event times: " << endl;
   //for ( std::vector<double>::iterator Pt_itr = Pt.begin();
   //      Pt_itr != Pt.end();
   //      ++Pt_itr)
   //{
   //   cout << *Pt_itr << ", ";
   //}
   //cout << endl;

   out_text_file_p.close();
   out_text_file_w.close();
   out_text_file_w_withdrift.close();
   out_text_file_w_withdrift2.close();

   out_text_file_p_stats.close();
   out_text_file_w_stats.close();
   out_text_file_w_withdrift_stats.close();
   out_text_file_w_withdrift2_stats.close();
   for (size_t ii= 0; ii < Ww.size(); ++ii) 
   {
      delete (Ww.back());
      Ww.pop_back();
   }
   for (size_t ii= 0; ii < Ww_withdrift.size(); ++ii) 
   {
      delete (Ww_withdrift.back());
      Ww_withdrift.pop_back();
   }
   for (size_t ii= 0; ii < Ww_withdrift2.size(); ++ii) 
   {
      delete (Ww_withdrift2.back());
      Ww_withdrift2.pop_back();
   }
   for (size_t ii= 0; ii < Pp.size(); ++ii) 
   {
      delete (Pp.back());
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
