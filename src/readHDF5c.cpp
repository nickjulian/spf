/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: readHDF5c.cpp
// Purpose:

#ifndef READHDF5C_CPP
#define READHDF5C_CPP

#include "readHDF5c.hpp"


int SPF_NS::read_phi_from_hdf5( 
      const hid_t inFile_id,
      std::vector<double>& phi, 
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
   hid_t phi_dataset_id, phi_dataspace_id;
   herr_t status; status = 0;
   hssize_t N_total;
   // open the hdf5 file
   phi_dataset_id = H5Dopen( inFile_id, datasetPath.c_str(), H5P_DEFAULT);
   phi_dataspace_id = H5Dget_space( phi_dataset_id );

   // reacquire spatial dimensions from the input data
   int ndims;
   ndims = H5Sget_simple_extent_ndims( phi_dataspace_id );
   //if (ndims != phi.size())
   //{
   //   //failflag += -1;
   //   if ( mynode == rootnode )
   //   {
   //      cout << "Error, ndims changed between file accesses ..." << endl;
   //   }
   //   H5Sclose( phi_dataspace_id );
   //   H5Dclose( phi_dataset_id );
   //   return EXIT_FAILURE;
   //}
   hsize_t dims[ndims];
   H5Sget_simple_extent_dims( phi_dataspace_id, dims, NULL);

   //cout << "node " << mynode << " initial phi values: "; // debug
   //for (std::vector<double>::const_iterator itr=phi.begin(); // debug
   //      itr != phi.end(); ++itr) // debug
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
      phi.resize( idx_end[0] - idx_start[0] +3 );
   else if ( ndims == 2 ) 
      phi.resize( (idx_end[0] - idx_start[0] +3) * dims[1] );
   else if ( ndims == 3 ) 
      phi.resize( (idx_end[0] - idx_start[0] +3) * dims[1] * dims[2]);
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

   //cout << "node " << mynode << " phi.size(): " << phi.size() << endl; // debug
   //hid_t phi_memspace_subset_id;
   //phi_memspace_subset_id = H5Screate_simple( ndims, dims, NULL);
   status = H5Sselect_hyperslab( 
                                 phi_dataspace_id,
                                 H5S_SELECT_SET, // H5S_seloper_t op
                                 offset, // const hsize_t *start,
                                 stride, // const hsize_t *stride,
                                 count, // const hsize_t *count,
                                 block // const hsize_t *block,
                                 );
   if ( status < 0 )
   {
      cout << "error: node " << mynode 
         << " failed to read phi from input file "//<< inputFileName
         << endl; // assumed to be executed only on one node
      H5Sclose( phi_dataspace_id );
      H5Dclose( phi_dataset_id );
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
         phi_dataset_id, // dataset_id,
         H5T_NATIVE_DOUBLE, // mem_type_id
         H5S_ALL, // mem_space_id
         phi_dataspace_id, // file_space_id
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

   // copy from the read buffer to phi
   for ( size_t ii=0; ii < idx_end[0] - idx_start[0] + 1; ++ii)
   {
      
      if ( ndims == 1)
         phi[ 1 + ii ] = input_buffer[idx_start[0] + ii];
      if ( ndims == 2)
      {
         for ( size_t jj=0; jj < dims[1]; ++jj)
            phi[ jj + dims[1]*(dims[1] + ii)] 
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
            //<< ", to phi[" <<  // debug
            //kk + dims[2]*(jj + dims[1]*(1 + ii))  // debug
            //<< "]" << endl; // debug

               phi[ 
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
         << " failed to read phi from input file "//<< inputFileName
         << endl; 
      H5Sclose( phi_dataspace_id );
      H5Dclose( phi_dataset_id );
      return EXIT_FAILURE;
   }

   H5Sclose( phi_dataspace_id );
   H5Dclose( phi_dataset_id );

   //cout << "node " << mynode << " phi values after reading: "; // debug
   //for (std::vector<double>::const_iterator itr=phi.begin(); // debug
   //      itr != phi.end(); ++itr) // debug
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
   hid_t phi_dataset_id, phi_dataspace_id;
   phi_dataset_id = H5Dopen2( inFile_id, datasetPath.c_str(), H5P_DEFAULT);
   if ( phi_dataset_id < 0 ) flags.fail = -1;
   if ( check_for_failure( flags, comm ) == true )
   {
      cout << "node " << mynode // debug 
         << " could not read '/phi' dataset" // debug
         << endl; // debug
      H5Dclose( phi_dataset_id );
      return EXIT_FAILURE;
   }
   phi_dataspace_id = H5Dget_space( phi_dataset_id );
   if ( phi_dataspace_id < 0 ) flags.fail = -1;
   if ( check_for_failure( flags, comm ) != 0)
   {
      cout << "node " << mynode  // debug
         << " could not read '/phi' dataspace" // debug
         << endl; // debug
      H5Sclose( phi_dataspace_id );
      H5Dclose( phi_dataset_id );
      return EXIT_FAILURE;
   }

   int ndims, ndims2;
   ndims = H5Sget_simple_extent_ndims( phi_dataspace_id );
   dims.resize(3, 1);
   hsize_t h5dims[ndims];
   ndims2 = H5Sget_simple_extent_dims( phi_dataspace_id, h5dims, NULL);
   if ( ndims2 < 1 || ndims2 > 3 ) flags.fail = -1;
   if ( check_for_failure( flags, comm ) == true )
   {
      cout << "node " << mynode  // debug
         << " could not read dimensions of phi_dataspace" // debug
         << endl; // debug
      H5Sclose( phi_dataspace_id );
      H5Dclose( phi_dataset_id );
      return EXIT_FAILURE;
   }
   
   //if (mynode == rootnode ) cout << "dimensions of phi: "; // debug
   for (size_t ii=0; ii < ndims2; ++ii)
   {
      dims[ii] = h5dims[ii];
      //if (mynode == rootnode ) cout << dims[ii] << ", "; // debug
   }
   if (mynode == rootnode) cout << endl; // debug

   H5Sclose( phi_dataspace_id );
   H5Dclose( phi_dataset_id );
   return EXIT_SUCCESS;
}

int SPF_NS::read_phi_from_hdf5_singlenode( 
      const hid_t inFile_id,
      const string& group_name,
      int& Nx, int& Ny, int& Nz,
      std::vector<double>& phi
      //const std::vector<int>& periodicity,
      )
{
   hid_t phi_dataset_id, phi_dataspace_id, group_id;
   herr_t status; status = 0;
   hssize_t N_total;
   // open the hdf5 file
   group_id = H5Gopen2(inFile_id, group_name.c_str(), H5P_DEFAULT);
   phi_dataset_id = H5Dopen2(group_id, "phi", H5P_DEFAULT);
   phi_dataspace_id = H5Dget_space( phi_dataset_id );

   // reacquire spatial dimensions from the input data
   int ndims;
   ndims = H5Sget_simple_extent_ndims( phi_dataspace_id );

   hsize_t dims[ndims];
   H5Sget_simple_extent_dims( phi_dataspace_id, dims, NULL);

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
   if ( phi.size() != total_number_of_elements )
   {
      cout << "resizing phi since phi.size() " << phi.size();
      if (ndims == 3)
      {
         Nx = dims[0];
         Ny = dims[1];
         Nz = dims[2];
         cout << phi.size() << " != dims[0]*dims[1]*dims[2] ";
      } 
      else if (ndims == 2) 
      { 
         cout << phi.size() << " != dims[0]*dims[1] ";
         Nx = dims[0];
         Ny = dims[1];
         Nz = 1;
      }
      else if (ndims == 1)
      {
         cout << phi.size() << " != dims[0] ";
         Nx = dims[0];
         Ny = 1;
         Nz = 1;
      }
      else
      {
         cout << "Error: ndims not among {1,2,3}" << endl;
         H5Sclose( phi_dataspace_id );
         H5Dclose( phi_dataset_id );
         H5Gclose( group_id );
         return EXIT_FAILURE;
      }
      cout << total_number_of_elements << endl;
   }

   phi.resize( total_number_of_elements );

   //std::vector<double> input_buffer(total_number_of_elements,0);
   status = H5Dread(
         phi_dataset_id, // dataset_id,
         H5T_NATIVE_DOUBLE, // mem_type_id
         H5S_ALL, // mem_space_id
         phi_dataspace_id, // file_space_id
         H5P_DEFAULT, // xfer_plist_id
         &phi[0] // void* buf
         //&input_buffer[0] // void* buf
         );

   if ( status < 0 )
   {
      cout << "error: " 
         << " failed to read phi from input file "//<< inputFileName
         << endl; 
      H5Sclose( phi_dataspace_id );
      H5Dclose( phi_dataset_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   H5Sclose( phi_dataspace_id );
   H5Dclose( phi_dataset_id );
   H5Gclose( group_id );

   return EXIT_SUCCESS;
}
#endif
