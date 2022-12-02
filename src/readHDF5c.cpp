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
// File: readHDF5c.cpp

#ifndef READHDF5C_CPP
#define READHDF5C_CPP

#include "readHDF5c.hpp"


int SPF_NS::read_dataset_from_hdf5(
      const hid_t inFile_id,
      std::vector<double>& data,
      const std::string datasetPath,
      const std::vector<size_t>& idx_start,
      const std::vector<size_t>& idx_end,
      //const std::vector<int>& periodicity,
      int_flags& flags,
      const int& mynode,
      const int& rootnode,
      const int& totalnodes,
      MPI_Comm comm
      )
{
   hid_t dataset_id, dataspace_id;
   herr_t status; status = 0;
   hssize_t N_total;
   // open the hdf5 file
   dataset_id = H5Dopen( inFile_id, datasetPath.c_str(), H5P_DEFAULT);
   dataspace_id = H5Dget_space( dataset_id );

   // reacquire spatial dimensions from the input data
   int ndims;
   ndims = H5Sget_simple_extent_ndims( dataspace_id );
   //if (ndims != data.size())
   //{
   //   //failflag += -1;
   //   if ( mynode == rootnode )
   //   {
   //      cout << "Error, ndims changed between file accesses ..." << endl;
   //   }
   //   H5Sclose( dataspace_id );
   //   H5Dclose( dataset_id );
   //   return EXIT_FAILURE;
   //}
   hsize_t dims[ndims];
   H5Sget_simple_extent_dims( dataspace_id, dims, NULL);

   //cout << "node " << mynode << " initial " << datasetPath  // debug
   // << " values: "; // debug
   //for (std::vector<double>::const_iterator itr=data.begin(); // debug
   //      itr != data.end(); ++itr) // debug
   //{ // debug
   //   cout << *itr << ", "; // debug
   //} // debug
   //cout << endl; // debug

   hsize_t offset[ndims];
   for ( size_t ii=0; ii < ndims; ++ii ) 
   {
      offset[ii] = idx_start[ii];
      //cout << "node " << mynode << " offset[" << ii << "]: " << offset[ii] << endl; // debug
   }
   hsize_t stride[ndims];
   for ( size_t ii=0; ii < ndims; ++ii ) stride[ii] = 1;
   hsize_t count[ndims];
   for ( size_t ii=0; ii < ndims; ++ii ) 
      count[ii] = idx_end[ii] - idx_start[ii] + 1;
   //for ( size_t ii=0; ii < ndims; ++ii ) 
   //{
   //   //count[ii] = () + 1;
   //   //cout << "node " << mynode << " count[" << ii << "]: " << count[ii] << endl; // debug
   //}
   hsize_t block[ndims];
   for ( size_t ii=0; ii < ndims; ++ii) block[ii] = 1;
   
   // resize the local data container to fit the input
   if ( ndims == 1 ) 
      data.resize( idx_end[0] - idx_start[0] +3 );
   else if ( ndims == 2 ) 
      data.resize( (idx_end[0] - idx_start[0] +3) * dims[1] );
   else if ( ndims == 3 ) 
      data.resize( (idx_end[0] - idx_start[0] +3) * dims[1] * dims[2]);
   else
   {
      if ( mynode == rootnode )
      {
         cout << "Error, spatial dimensions of input data out of range: "
            << "ndims == " << ndims
            << endl;
      }
      return EXIT_FAILURE;
   }

   //cout << "node " << mynode << " data.size(): " << data.size() // debug
   // << endl; // debug
   //hid_t data_memspace_subset_id;
   //data_memspace_subset_id = H5Screate_simple( ndims, dims, NULL);
   status = H5Sselect_hyperslab( 
                                 dataspace_id,
                                 H5S_SELECT_SET, // H5S_seloper_t op
                                 offset, // const hsize_t *start,
                                 stride, // const hsize_t *stride,
                                 count, // const hsize_t *count,
                                 block // const hsize_t *block,
                                 );
   if ( status < 0 )
   {
      cout << "error: node " << mynode 
         << " failed to read " << datasetPath << " from input file "
         //<< inputFileName
         << endl; // TODO: send and write the error from only root node
      H5Sclose( dataspace_id );
      H5Dclose( dataset_id );
      return EXIT_FAILURE;
   }

   // It seems that the hdf5 interface requires the destination of 
   //  a read operation to be the size of the entire dataspace, 
   //  not just the subset selected by the H5Sselect_hyperslab command.
   int total_number_of_elements;
   if ( ndims == 1 ) total_number_of_elements = dims[0];
   else if ( ndims == 2 ) total_number_of_elements = dims[0] * dims[1];
   else if ( ndims == 3 ) 
      total_number_of_elements = dims[0] * dims[1] * dims[2];

   std::vector<double> input_buffer(total_number_of_elements,0);
   status = H5Dread(
         dataset_id, // dataset_id,
         H5T_NATIVE_DOUBLE, // mem_type_id
         H5S_ALL, // mem_space_id
         dataspace_id, // file_space_id
         H5P_DEFAULT, // xfer_plist_id
         &input_buffer[0] // void* buf
         );
   //cout << "node " << mynode << " input_buffer: "; // debug
   //for (std::vector<double>::const_iterator itr=input_buffer.begin();//debug
   //      itr != input_buffer.end(); ++itr) // debug
   //{ // debug
   //   cout << *itr << ", "; // debug
   //} // debug
   //cout << endl; // debug

   // copy from the read buffer to data variable
   for ( size_t ii=0; ii < idx_end[0] - idx_start[0] + 1; ++ii)
   {
      
      if ( ndims == 1)
         data[ 1 + ii ] = input_buffer[idx_start[0] + ii];
      if ( ndims == 2)
      {
         for ( size_t jj=0; jj < dims[1]; ++jj)
            data[ jj + dims[1]*(dims[1] + ii)] 
               = input_buffer[jj + dims[1]*(idx_start[0] + ii)];
      }
      if ( ndims == 3)
      {
         for ( size_t jj=0; jj < dims[1]; ++jj)
            for ( size_t kk=0; kk < dims[2]; ++kk)
            {
            //cout << "assigning input_buffer[" <<  // debug
            //kk + dims[2]*(jj + dims[1]*(idx_start[0] + ii))
            //<< "]: " <<  // debug
            //input_buffer[
            //kk + dims[2]*(jj + dims[1]*(idx_start[0] + ii))
            //] // debug
            //<< ", to data[" <<  // debug
            //kk + dims[2]*(jj + dims[1]*(1 + ii))  // debug
            //<< "]" << endl; // debug

               data[ 
                  kk + dims[2]*(jj + dims[1]*(1 + ii)) 
               ]
                  = input_buffer[
                     kk + dims[2]*(jj + dims[1]*(idx_start[0] + ii))
                     ];
            }
      }
   }

   if ( status < 0 )
   {
      cout << "error: node " << mynode 
         << " failed to read " << datasetPath << " from input file "
         //<< inputFileName
         << endl; 
      H5Sclose( dataspace_id );
      H5Dclose( dataset_id );
      return EXIT_FAILURE;
   }

   H5Sclose( dataspace_id );
   H5Dclose( dataset_id );

   //cout << "node " << mynode << " data values after reading: "; // debug
   //for (std::vector<double>::const_iterator itr=data.begin(); // debug
   //      itr != data.end(); ++itr) // debug
   //{ // debug
   //   cout << *itr << ", "; // debug
   //} // debug
   //cout << endl; // debug

   return EXIT_SUCCESS;
}

int SPF_NS::determine_local_idxs(
      const std::vector<hsize_t>& dims,
      const int& mynode,
      const int& rootnode,
      const int& totalnodes,
      int& Nx_local,
      std::vector<size_t>& idx_start,
      std::vector<size_t>& idx_end
      )
{
   int ndims; ndims = dims.size();
   if (ndims <= 0)
   {
      return EXIT_FAILURE;
   }

   int Nx; Nx = dims[0];
   if ( totalnodes > Nx )
   {
      if ( mynode == rootnode )
      {
         cout << "Error: number of compute nodes"
               << " > # of points on phi domain ..." << endl;
      }
      return EXIT_FAILURE;
   }

   // This function assigns indices of non-ghost local data to use in a
   //  global data array.
   if ( Nx % totalnodes == 0 )
   {
      Nx_local = Nx / totalnodes; // count of non-ghost rows along x
      idx_start[0] = mynode * Nx_local ; // idx of non-ghost data beginning
      idx_end[0] = (mynode+1) * Nx_local -1;
   }
   else
   {
      int remainderNx; remainderNx = Nx % totalnodes;
      int commonNx; commonNx = (Nx - remainderNx) / totalnodes;
      if ( mynode < remainderNx )
      {
         Nx_local = 1 + commonNx;
         idx_start[0] = mynode * Nx_local;
         idx_end[0] = (mynode +1)* Nx_local -1;
      }
      else
      {
         Nx_local = commonNx;
         idx_start[0] = remainderNx * (Nx_local +1)
                        + (mynode - (remainderNx -1) -1) * Nx_local;
         idx_end[0] = idx_start[0] + Nx_local -1;
      }
   }
   for ( size_t ii=1; ii < ndims; ++ii)
   {
      idx_start[ii] = 0;
      idx_end[ii] = dims[ii] - 1;
   }

   return EXIT_SUCCESS;
}

int SPF_NS::read_dims_from_hdf5( 
         const hid_t inFile_id,
         const std::string datasetPath,
         std::vector<hsize_t>& dims,
         int_flags& flags,
         const int& mynode,
         const int& rootnode,
         const int& totalnodes,
         MPI_Comm comm
      )
{
   // resizes dims and sets dims[i] = # elements in i^{th} dimension
   // assumes all field data have same dimensions
   hid_t dataset_id, dataspace_id;
   dataset_id = H5Dopen2( inFile_id, datasetPath.c_str(), H5P_DEFAULT);
   if ( dataset_id < 0 ) flags.fail = -1;
   if ( check_for_failure( flags, comm ) == true )
   {
      cout << "node " << mynode // debug 
         << " could not read " << datasetPath << " dataset" // debug
         << endl; // debug
      H5Dclose( dataset_id );
      return EXIT_FAILURE;
   }
   dataspace_id = H5Dget_space( dataset_id );
   if ( dataspace_id < 0 ) flags.fail = -1;
   if ( check_for_failure( flags, comm ) != 0)
   {
      cout << "node " << mynode  // debug
         << " could not read " << datasetPath << " dataspace" // debug
         << endl; // debug
      H5Sclose( dataspace_id );
      H5Dclose( dataset_id );
      return EXIT_FAILURE;
   }

   int ndims, ndims2;
   ndims = H5Sget_simple_extent_ndims( dataspace_id );
   dims.resize(3, 1);
   hsize_t h5dims[ndims];
   ndims2 = H5Sget_simple_extent_dims( dataspace_id, h5dims, NULL);
   if ( ndims2 < 1 || ndims2 > 3 ) flags.fail = -1;
   if ( check_for_failure( flags, comm ) == true )
   {
      cout << "node " << mynode  // debug
         << " could not read dimensions of " << datasetPath // debug
         << " dataspace" // debug
         << endl; // debug
      H5Sclose( dataspace_id );
      H5Dclose( dataset_id );
      return EXIT_FAILURE;
   }
   
   //if (mynode == rootnode ) cout << "dimensions of " // debug 
   // << datasetPath << ": "; // debug
   for (size_t ii=0; ii < ndims2; ++ii)
   {
      dims[ii] = h5dims[ii];
      //if (mynode == rootnode ) cout << dims[ii] << ", "; // debug
   }
   if (mynode == rootnode) cout << endl; // debug

   H5Sclose( dataspace_id );
   H5Dclose( dataset_id );
   return EXIT_SUCCESS;
}

int SPF_NS::read_dataset_from_hdf5_singlenode( 
      const hid_t inFile_id,
      const string& group_name,
      const string& datasetPath,
      int& Nx, int& Ny, int& Nz,
      std::vector<double>& data
      //const std::vector<int>& periodicity,
      )
{
   hid_t dataset_id, dataspace_id, group_id;
   herr_t status; status = 0;
   hssize_t N_total;
   // open the hdf5 file
   group_id = H5Gopen2(inFile_id, group_name.c_str(), H5P_DEFAULT);
   dataset_id = H5Dopen2(group_id, datasetPath.c_str(), H5P_DEFAULT);
   dataspace_id = H5Dget_space( dataset_id );

   // reacquire spatial dimensions from the input data
   int ndims;
   ndims = H5Sget_simple_extent_ndims( dataspace_id );

   hsize_t dims[ndims];
   H5Sget_simple_extent_dims( dataspace_id, dims, NULL);

   hsize_t stride[ndims];
   for ( size_t ii=0; ii < ndims; ++ii ) stride[ii] = 1;

   hsize_t count[ndims];
   for ( size_t ii=0; ii < ndims; ++ii ) 
      count[ii] = dims[ii];
      //count[ii] = idx_end[ii] - idx_start[ii] + 1;

   hsize_t block[ndims];
   for ( size_t ii=0; ii < ndims; ++ii) block[ii] = 1;
   
   int total_number_of_elements;
   if ( ndims == 1 ) total_number_of_elements = dims[0];
   else if ( ndims == 2 ) total_number_of_elements = dims[0] * dims[1];
   else if ( ndims == 3 ) 
      total_number_of_elements = dims[0] * dims[1] * dims[2];

   // resize the local data container to fit the input
   if ( data.size() != total_number_of_elements )
   {
      cout << "resizing data since data.size() " << data.size();
      if (ndims == 3)
      {
         Nx = dims[0];
         Ny = dims[1];
         Nz = dims[2];
         cout << data.size() << " != dims[0]*dims[1]*dims[2] ";
      } 
      else if (ndims == 2) 
      { 
         cout << data.size() << " != dims[0]*dims[1] ";
         Nx = dims[0];
         Ny = dims[1];
         Nz = 1;
      }
      else if (ndims == 1)
      {
         cout << data.size() << " != dims[0] ";
         Nx = dims[0];
         Ny = 1;
         Nz = 1;
      }
      else
      {
         cout << "Error: ndims not among {1,2,3}" << endl;
         H5Sclose( dataspace_id );
         H5Dclose( dataset_id );
         H5Gclose( group_id );
         return EXIT_FAILURE;
      }
      cout << total_number_of_elements << endl;
   }

   data.resize( total_number_of_elements );

   //std::vector<double> input_buffer(total_number_of_elements,0);
   status = H5Dread(
         dataset_id, // dataset_id,
         H5T_NATIVE_DOUBLE, // mem_type_id
         H5S_ALL, // mem_space_id
         dataspace_id, // file_space_id
         H5P_DEFAULT, // xfer_plist_id
         &data[0] // void* buf
         //&input_buffer[0] // void* buf
         );

   if ( status < 0 )
   {
      cout << "error: " 
         << " failed to read " << datasetPath << "from input file "
         //<< inputFileName
         << endl; 
      H5Sclose( dataspace_id );
      H5Dclose( dataset_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   H5Sclose( dataspace_id );
   H5Dclose( dataset_id );
   H5Gclose( group_id );

   return EXIT_SUCCESS;
}
#endif
