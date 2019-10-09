/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: readHDF5c.hpp
// Purpose:

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
   hid_t outFileAverage_id, outFileVariance_id;
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // open files to write the statistical evolution to
   string outputFileNameVariance; 
   outputFileNameVariance = output_prefix + "_variance.h5";
   outFileVariance_id = H5Fcreate(outputFileNameVariance.c_str(), 
                        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   string outputFileNameAverage; 
   outputFileNameAverage = output_prefix + "_average.h5";
   outFileAverage_id = H5Fcreate(outputFileNameAverage.c_str(), 
                        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   if ((outFileVariance_id < 0) || (outFileAverage_id < 0)) 
   {
         cout << "Error, failed to open files for writing: " 
            << endl << "   "
            << outputFileNameVariance << endl
            << endl << "   "
            << outputFileNameAverage << endl;
      H5Fclose( outFileVariance_id );
      H5Fclose( outFileAverage_id );
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
         //group_name = "Step#" + std::to_string(time_step);
         // convert time_step from into to string
         std::ostringstream sstime_step;
         sstime_step << time_step;
         group_name = "Step#" + sstime_step.str();
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
         cout << "calculating average and variance at each voxel" // debug
            << endl; // debug
         for ( size_t ii=0; ii < Nx; ++ii)
            for ( size_t jj=0; jj < Ny; ++jj)
               for ( size_t kk=0; kk < Nz; ++kk)
               {
                  phi_sum[kk + Nz*(jj + Ny*ii)] 
                   = (phi_sum[kk + Nz*(jj + Ny*ii)] )/ Nstates;

                  phi_sqr_sum[kk + Nz*(jj + Ny*ii)] 
                     = phi_sqr_sum[kk + Nz*(jj + Ny*ii)] / Nstates
                        - (phi_sum[kk + Nz*(jj + Ny*ii)]
                           * phi_sum[kk + Nz*(jj + Ny*ii)]);
               }
         cout << "saving averages to file "  // debug
            << outputFileNameAverage << endl; // debug
         // write the average and variance of the current time to file
         append_phi_to_hdf5_singlenode( 
               outFileAverage_id,
               time_step,
               time,
               "average",
               dims, // dims.size() == 3
               phi_sum
               );

         cout << "saving variances to file "  // debug
            << outputFileNameVariance << endl; // debug
         append_phi_to_hdf5_singlenode( 
               outFileVariance_id,
               time_step,
               time,
               "variance",
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

   H5Fclose( outFileVariance_id );
   H5Fclose( outFileAverage_id );
   MPI_Finalize();
   return EXIT_SUCCESS;
}

#endif
