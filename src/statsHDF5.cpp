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
// File: readHDF5c.hpp

#ifndef STATSHDF5C_CPP
#define STATSHDF5C_CPP

#include <iostream>  // cout, cin, cerr, endl
#include <iomanip>   // setw, setprecision
#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE
#include <string>
#include <sstream>   // ostringstream
#include <math.h> // floor
#include <vector>
#include <mpi.h>
#include "../include/hdf5.h"  

#include "read_cmdline.hpp"
#include "read_parameter_files.hpp"
#include "readHDF5c.hpp"
#include "writeHDF5c.hpp"

using std::cout;
using std::endl;
using std::string;
using std::cerr;

using namespace SPF_NS;

int main( int argc, char* argv[])
{
   int failflag; failflag = 0;

   std::vector<hsize_t> dims;
   int Nx, Ny, Nz, Nt;
   Nx = 1; Ny = 1; Nz = 1;
   int time_step; time_step = 0;
   int write_period; write_period = 1;
   double time, dt; time = 0; dt = 1;
   double rate_scale_factor; rate_scale_factor = 1.0;

   std::string inputFileName;
   std::vector<string> args( argv, argv + argc );
   std::string output_prefix;
   bool flag_calcstat;
   int mynode;
   int rootnode; rootnode = 0;
   int totalnodes;
   MPI_Init( &argc, &argv);
   MPI_Comm_size( MPI_COMM_WORLD, &totalnodes );
   MPI_Comm_rank( MPI_COMM_WORLD, &mynode );
   if ( mynode != rootnode ) 
   {
      MPI_Finalize();
      return EXIT_SUCCESS;
   }
   if ( read_cmdline_options(
            args,
            dt,
            Nt,
            rate_scale_factor,
            write_period,
            flag_calcstat,
            output_prefix,
            inputFileName,
            mynode,
            rootnode,
            MPI_COMM_WORLD
            ) != EXIT_SUCCESS )
   {
      if ( mynode == rootnode ) PRINT_USAGE;
      MPI_Finalize();
      return EXIT_FAILURE;
   }

   std::vector<std::string> state_file_names;

   if ( read_state_file_list(
            inputFileName,
            state_file_names
            ) != EXIT_SUCCESS)
   {
      cerr << "Error, failed to read " << inputFileName << endl;
      MPI_Finalize();
      return EXIT_FAILURE;
   }
   size_t Nstates; Nstates = state_file_names.size();

   ////////////////////////////////////////////////////////////////////
   // variables reused over time steps and state files
   string group_name;
   std::vector<double> phi;
   std::vector<double> phi_sum;
   std::vector<double> phi_sqr_sum;
   hid_t inFile_id, phi_dataset_id, phi_dataspace_id, phi_group_id;
   hid_t outFileSum_id, outFileSqrSum_id;
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // open files to write the statistical evolution to
   string outputFileNameSqrSum; 
   outputFileNameSqrSum = output_prefix + "_sqr_sum.h5";
   outFileSqrSum_id = H5Fcreate(outputFileNameSqrSum.c_str(), 
                        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   string outputFileNameSum; 
   outputFileNameSum = output_prefix + "_sum.h5";
   outFileSum_id = H5Fcreate(outputFileNameSum.c_str(), 
                        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   if ((outFileSqrSum_id < 0) || (outFileSum_id < 0)) 
   {
         cout << "Error, failed to open files for writing: " 
            << endl << "   "
            << outputFileNameSqrSum << endl
            << endl << "   "
            << outputFileNameSum << endl;
      H5Fclose( outFileSqrSum_id );
      H5Fclose( outFileSum_id );
      return EXIT_FAILURE;
   }
   ////////////////////////////////////////////////////////////////////


   /*--------------------------------------------------------------*/
   /* begin loop over time ----------------------------------------*/
   /*--------------------------------------------------------------*/
   for (time_step = 0; time_step <= Nt; ++time_step)
   {
      if ( time_step % write_period == 0 )
      {
         //group_name = "Step" + std::to_string(time_step);
         // convert time_step from into to string
         std::ostringstream sstime_step;
         sstime_step << time_step;
         group_name = "Step" + sstime_step.str();
         cout << "reading group " << group_name << endl; // debug

         /*-----------------------------------------------------------*/
         /* begin loop over state files ------------------------------*/
         /*-----------------------------------------------------------*/
         for (std::vector<string>::const_iterator 
                  itr = state_file_names.begin(); 
                  itr != state_file_names.end(); ++itr) 
         {
            // assign inputFileName
            inputFileName = *itr;

            //////////////////////////////////////////////////////////////
            // open input data file
            inFile_id = H5Fopen(inputFileName.c_str(), 
                              H5F_ACC_RDONLY, H5P_DEFAULT);

            if (inFile_id < 0) 
            {
               cout << "Error, failed to open input file: " 
                  << inputFileName << endl;
               H5Fclose( inFile_id );
               return EXIT_FAILURE;
            }
            //////////////////////////////////////////////////////////////

            if ( read_phi_from_hdf5_singlenode( 
                     inFile_id,
                     group_name,
                     Nx, Ny, Nz,
                     phi
                     ) == EXIT_FAILURE)
            {
               cout << "Error, could not read /" << group_name 
                  << "/phi from file " << inputFileName << endl;
               H5Fclose( inFile_id );
               MPI_Finalize();
               return EXIT_FAILURE;
            }
            H5Fclose( inFile_id );

            dims.resize(3);
            dims[0] = Nx;
            dims[1] = Ny;
            dims[2] = Nz;

            // this resize should only happen the first loop
            if ( phi_sum.size() != Nx * Ny * Nz )
            {
               cout << "resizing phi_sum to match phi" << endl;
               phi_sum.resize(Nx * Ny * Nz, 0);
            }
            if ( phi_sqr_sum.size() != Nx * Ny * Nz )
            {
               cout << "resizing phi_sqr_sum to match phi" << endl;
               phi_sqr_sum.resize(Nx * Ny * Nz, 0);
            }

            // loop over voxels
            cout << "looping over voxels at time step " // debug
               << time_step << endl; // debug
            for ( size_t ii=0; ii < Nx; ++ii)
               for ( size_t jj=0; jj < Ny; ++jj)
                  for ( size_t kk=0; kk < Nz; ++kk)
                  {
                     // accumulate the average and average of squares
                     phi_sum[kk + Nz*(jj + Ny*ii)] 
                        += phi[kk + Nz*(jj + Ny*ii)];

                     phi_sqr_sum[kk + Nz*(jj + Ny*ii)] 
                        += phi[kk + Nz*(jj + Ny*ii)] 
                           * phi[kk + Nz*(jj + Ny*ii)];
                  }
            cout << "done acumulating sums for this time step" // debug
               << endl; // debug
         }
         /*--------------------------------------------------------*/
         /* end loop over state files -----------------------------*/
         /*--------------------------------------------------------*/

         // calculate the average and variance of the current time step
         //cout << "calculating average and variance at each voxel" // debug
         //   << endl; // debug
         //for ( size_t ii=0; ii < Nx; ++ii)
         //   for ( size_t jj=0; jj < Ny; ++jj)
         //      for ( size_t kk=0; kk < Nz; ++kk)
         //      {
         //         phi_sum[kk + Nz*(jj + Ny*ii)] 
         //          = (phi_sum[kk + Nz*(jj + Ny*ii)] )/ Nstates;

         //         phi_sqr_sum[kk + Nz*(jj + Ny*ii)] 
         //            = phi_sqr_sum[kk + Nz*(jj + Ny*ii)] / Nstates
         //               - (phi_sum[kk + Nz*(jj + Ny*ii)]
         //                  * phi_sum[kk + Nz*(jj + Ny*ii)]);
         //      }
         cout << "saving sums to file "  // debug
            << outputFileNameSum << endl; // debug
         // write the average and variance of the current time to file
         append_phi_to_hdf5_singlenode( 
               outFileSum_id,
               time_step,
               time,
               "sum",
               dims, // dims.size() == 3
               phi_sum
               );

         cout << "saving sum of squares to file "  // debug
            << outputFileNameSqrSum << endl; // debug
         append_phi_to_hdf5_singlenode( 
               outFileSqrSum_id,
               time_step,
               time,
               "sqr_sum",
               dims, // dims.size() == 3
               phi_sqr_sum
               );

         // reset sums
         for ( size_t ii=0; ii < Nx; ++ii)
            for ( size_t jj=0; jj < Ny; ++jj)
               for ( size_t kk=0; kk < Nz; ++kk)
               {
                  phi_sum[kk + Nz*(jj + Ny*ii)] = 0.0;
                  phi_sqr_sum[kk + Nz*(jj + Ny*ii)] = 0.0;
               }
      } // if ( time_step % write_period == 0)
   }
   /*--------------------------------------------------------------*/
   /* end loop over time ------------------------------------------*/
   /*--------------------------------------------------------------*/

   H5Fclose( outFileSqrSum_id );
   H5Fclose( outFileSum_id );
   MPI_Finalize();
   return EXIT_SUCCESS;
}

#endif
