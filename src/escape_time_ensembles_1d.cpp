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
#include "flags.hpp"

// TODO: execute an ensemble of trajectories having uniformly distributed initial positions, ending each path when it reaches some maximum value and save these first exit times to a file for external analysis

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
   int_flags flags;
   string output_prefix;
   
   // increments
   double dt, sqrtdt; 
   size_t Nensemble; Nensemble = 100000;
   //size_t Nensemble; Nensemble = 1000;
   //size_t Nensemble; Nensemble = 10;
   //size_t Nensemble; Nensemble = 1;
   int Nt;
   //double dx, dy, dz;

   std::vector<string> args( argv, argv + argc );

   if (read_cmdline_options_for_escape_time_ensembles(
         args,
         dt,
         Nt,
         Nensemble,
         flags,
         output_prefix
         //mynode,
         //rootnode,
         //MPI_COMM_WORLD
         ) != EXIT_SUCCESS )
   {
      if ( mynode == rootnode ) PRINT_USAGE;

      MPI_Finalize();
      return EXIT_FAILURE;
   }

   double end_time; end_time = Nt * dt;
   cout << "end_time : " << end_time;
   double w_rate; w_rate = rate_scale_factor;
   double p_rate; p_rate = rate_scale_factor;
   //double p_rate_mu0; p_rate_mu0 = rate_scale_factor;
   double p_rate_mu0; p_rate_mu0 = 0.5* rate_scale_factor;
   double w_sqrtrate; 

   double drift; drift = w_rate;

   sqrtdt = sqrt(dt);
   w_sqrtrate = sqrt(w_rate);

   // random generator
   SPF_NS::random rr;

   // stochastic distributions
   //std::poisson_distribution<int> pd( p_rate );
   std::exponential_distribution<double> ed( p_rate );
   std::exponential_distribution<double> ed_mu0( p_rate_mu0 );
   std::bernoulli_distribution bd( 0.5 );

   // time variables
   std::vector<double> Pt;
   std::vector<double> Wt;
   std::vector<double> Pt_mu0;
   std::vector<double> Wt_drift;

   //std::list<std::vector<double>*> Pp;
   //std::list<std::vector<double>*> Xp; // dX(t) is the SDE driven by dP(t)
   //std::list<std::vector<double>*> Xp_exact;// path according to analysis
   
   // Poisson processes
   double p_dt;
   double p_f_dt, p_b_dt;
   int dP; dP = 1.0;
   int dP_f, dP_b;

   // Wiener processes
   double dW;
   std::normal_distribution<double> gd( 0.0, 1.0);//); // (mean, stddev)

   //std::list<std::vector<double>*> Ww; // W(t) is the stochastic process
   //std::list<std::vector<double>*> Xw; // dX(t) is the SDE driven by dW(t)
   //std::list<std::vector<double>*> Xw_exact;// path according to analysis

   //std::vector<double>::iterator Tt_itr;
   //std::vector<double>* Ww_ptr;
   //std::vector<double>* Pp_ptr;
   //std::vector<double> Xx_ptr;
   //std::vector<double> Xx_exact_ptr;


   //cout << " rate : "  << rate << endl; // debug

   double w_init; w_init = 5.0;
   double w_lb; w_lb = 0.0;
   double w_ub; w_ub = 10.0;

   int p_init; p_init = 0;
   int p_lb; p_lb = 0;
   int p_ub; p_ub = 10;

   double Ww; 
   double Ww_drift; 
   int Pp;
   int Pp_mu0;

   // loop over ensemble members
   for ( size_t en_count=0; en_count < Nensemble; ++en_count)
   {
      //Ww.push_back( new std::vector<double> );
      //Pp.push_back( new std::vector<double> );
      //Xx.push_back( new std::vector<double> );
      //Xx_exact.push_back( new std::vector<double> );

      //Ww_ptr = Ww.back();
      //Pp_ptr = Pp.back();
      //Xx_ptr = Xx.begin();
      //Xx_exact_ptr = Xx_exact.begin();

      // initial coordinate
      //Ww_ptr->push_back( 0.0 );
      //Pp_ptr->push_back( 0.0 );
      //Wt.push_back( 0.0 );
      //Pt.push_back( 0.0 );
      //Xx_ptr->push_back( 0.0 );
      //Xx_exact_ptr->push_back( 0.0 );

      // wiener process, loop over time
      Ww = w_init;
      for( int tt=0; tt < Nt; ++tt)
      {
         // TODO: evaluate the analytical solution
         // The SDE
         // dX(t) = a X(t) dt + \mu X(t) dW(t)
         // has exact solution
         // X(t) = X(0) exp( (a - 0.5 \mu^2)t + \mu W(t))
         

         // Wiener process
         // dW ~ \sqrt(\Delta t) * N(0, 1.0)
         // W_t = W_{t-1} + f(x) dt + g(x) dW + 0.5 g(x) g'(x) ((dW)^2 -dt)
         // g(x) = rate_scale_factor * ?something, which here is 1?
         // g'(x) = d/dx g(x) = 0
         // f(x) = drift
         dW = gd( rr.generator );
         Ww += (sqrtdt 
                 //* dW); // + 0, since g'(x) =0
                 * w_sqrtrate * dW); // + 0, since g'(x) =0

         // dW ~ \sqrt(\Delta t) * N( 0, 1.0) 

         if ( Ww > w_ub || Ww < w_lb ) 
         {
            Wt.push_back( tt );
            break;
         }
      } // loop over time

      Ww_drift = double(p_init);
      for( int tt=0; tt < Nt; ++tt)
      {
         // Wiener process with drift
         // dW ~ \sqrt(\Delta t) * N(0, 1.0)
         // W_t = W_{t-1} + f(x) dt + g(x) dW + 0.5 g(x) g'(x) ((dW)^2 -dt)
         // g(x) = rate_scale_factor * ?something, which here is 1?
         // g'(x) = d/dx g(x) = 0
         // f(x) = drift
         dW = gd( rr.generator );
         Ww_drift += (drift * dt) + (sqrtdt 
                 * w_sqrtrate * dW); // + 0, since g'(x) =0

         if ( Ww_drift > p_ub )//|| Ww_drift < p_lb )
         {
            Wt_drift.push_back( tt );
            break;
         }
      } // loop over time


      // Poisson process
      // loop over time until tt > end_time
      double tt; tt = 0.0;
      Pp = p_init;
      while ( tt <= end_time )
      {
         // Poisson process
         
         // sample time increment
         p_dt = ed( rr.generator );
         
         //if ( bd( rr.generator ) )
         //{
            dP = 1;
         //}
         //else
         //{
         //   dP = -1;
         //}
         Pp += dP;
         //Pp_ptr->push_back( (*Pp_ptr)[ Pp_ptr->size() -1]
         //                     + (p_dt * dP) );
         tt  += p_dt;
         if ( Pp > p_ub ) //|| Pp < p_lb ) 
         {
            Pt.push_back( tt );
            break;
         }
      }

      tt = 0.0;
      Pp_mu0 = int(std::round(w_init));
      while ( tt <= end_time )
      {
         // two Poisson processes in opposite directions
         
         // sample time increment
         p_f_dt = ed_mu0( rr.generator );
         p_b_dt = ed_mu0( rr.generator );
         
         if ( p_f_dt <= p_b_dt ) 
         {
            dP = 1.0;
            tt  += p_f_dt;
         }
         else
         {
            dP = -1.0;
            tt  += p_b_dt;
         }

         Pp_mu0 += dP;

         if ( Pp_mu0 > w_ub || Pp_mu0 < w_lb ) 
         {
            Pt_mu0.push_back( tt );
            break;
         }
      }

      //cout << "path # " << en_count << ", Pp_ptr->size() : " << Pp_ptr->size() << endl; // debug
   } // loop over ensemble members

   // text file output
   ofstream out_text_file_w;
   out_text_file_w.open(output_prefix + "_Nens" + std::to_string(Nensemble)
         + "_1d_escape_times_w.txt", ios::app | ios::ate);

   ofstream out_text_file_p;
   out_text_file_p.open(output_prefix + "_Nens" + std::to_string(Nensemble)
         + "_1d_escape_times_p.txt", ios::app | ios::ate);

   ofstream out_text_file_w_drift;
   out_text_file_w_drift.open(output_prefix + "_Nens" + std::to_string(Nensemble)
         + "_1d_escape_times_w_drift.txt", ios::app | ios::ate);

   ofstream out_text_file_p_mu0;
   out_text_file_p_mu0.open(output_prefix + "_Nens" + std::to_string(Nensemble)
         + "_1d_escape_times_p_mu0.txt", ios::app | ios::ate);

   if ( ! out_text_file_p.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix + "_Nens" + std::to_string(Nensemble) 
               +  "_1d_escape_times_p.txt" << endl;
      out_text_file_p.close();
      return EXIT_FAILURE;
   }

   if ( ! out_text_file_w.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix + "_Nens" + std::to_string(Nensemble) 
               +  "_1d_escape_times_w.txt" << endl;
      out_text_file_w.close();
      return EXIT_FAILURE;
   }

   if ( ! out_text_file_p_mu0.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix + "_Nens" + std::to_string(Nensemble) 
               +  "_1d_escape_times_p_mu0.txt" << endl;
      out_text_file_p_mu0.close();
      return EXIT_FAILURE;
   }

   if ( ! out_text_file_w_drift.good() )
   {
      cout << "warning: could not open output file: "
            << output_prefix + "_Nens" + std::to_string(Nensemble) 
               +  "_1d_escape_times_w_drift.txt" << endl;
      out_text_file_w_drift.close();
      return EXIT_FAILURE;
   }


   // write first passage times to text files
   for (std::vector<double>::iterator 
            itr = Wt.begin(); itr != Wt.end(); ++itr)
   {
      out_text_file_w << setw(12) << setprecision(5) << *itr << " ";
   }

   out_text_file_w << endl;

   for (std::vector<double>::iterator 
            itr = Pt.begin(); itr != Pt.end(); ++itr)
   {
      out_text_file_p << setw(12) << setprecision(5) << *itr << " ";
   }

   out_text_file_p << endl;

   for (std::vector<double>::iterator 
            itr = Wt_drift.begin(); itr != Wt_drift.end(); ++itr)
   {
      out_text_file_w_drift << setw(12) << setprecision(5) << *itr << " ";
   }

   out_text_file_w_drift << endl;

   for (std::vector<double>::iterator 
            itr = Pt_mu0.begin(); itr != Pt_mu0.end(); ++itr)
   {
      out_text_file_p_mu0 << setw(12) << setprecision(5) << *itr << " ";
   }

   out_text_file_p_mu0 << endl;

   out_text_file_p.close();
   out_text_file_w.close();
   out_text_file_p_mu0.close();
   out_text_file_w_drift.close();

   return EXIT_SUCCESS;
}
