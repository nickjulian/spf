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
   // read input parameters
   string outFilePrefix;
   std::vector<string> args( argv, argv + argc );
   if ( args.size() > 2 )
   {
      cout << "Usage: " << argv[0] << " <output file prefix>" << endl;
      return EXIT_FAILURE;
   }
   else if ( args.size() == 2 )
   {
      outFilePrefix = string( args[1] ) + "_";
   }

   int Npoints;
   Npoints  = 10000;

   // increments
   double dt, sqrtdt; dt = 1; sqrtdt = sqrt(dt);
   double dx, dy, dz;
   double rate1; rate1 = 2.0; // Poisson procces rate
   double rate2; rate2 = 1.0; // Poisson procces rate

   TEM_NS::random rr1( rate1, dt);
   TEM_NS::random rr2( rate2, dt);

   // Poisson processes
   std::vector<double> Pt;
   
   std::ostringstream Npoints_ostringstream;
   Npoints_ostringstream << Npoints;

   outFilePrefix = outFilePrefix + "exponential_times_" 
                     + Npoints_ostringstream.str() + "_time_increments.h5";

   int Nmisses; Nmisses = 0;
   //double poisson_rate1; poisson_rate1 = rr1.poisson_event_count.param();
   //double poisson_rate2; poisson_rate2 = rr2.poisson_event_count.param();
   int count1, count2; count1 = 0; count2 = 0;
   for( size_t p=0; p < Npoints; ++p)
   {
      //// Poisson trials, useful if at most one event occurs per interval
      //if ( rr.poisson_trial() )
      //{
      //   Pt.push_back( p * dt ); 
      //   Nmisses = 0;
      //}

      // sample a Poisson distribution to tally the number of events this dt
      //poisson_rate1 = 0.5 * poisson_rate1;
      //rr1.poisson_event_count.param( poisson_rate1 );
      //rr2.poisson_event_count.param( poisson_rate2 );
      count1 = rr1.poisson_event_count( rr1.generator ); 
      count2 = rr2.poisson_event_count( rr2.generator ); 
      Pt.push_back( count1 - count2 );
   }

   // output Pt to a file
   output_vector_to_hdf5( Pt, outFilePrefix);

   // clear the poisson times
   Pt.clear();

   return EXIT_SUCCESS;
}
