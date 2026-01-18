/* ----------------------------------------------------------------------
    SPF - Stochastic Phase Field
    Copyright (C) 2025
    Peng Geng <penggeng@g.ucla.edu>
    Nicholas Huebner Julian <njulian@ucla.edu>
    
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
// File: writeHDF5c.cpp

#ifndef WRITEHDF5C_CPP
#define WRITEHDF5C_CPP

#include "writeHDF5c.hpp"

int SPF_NS::append_fields_to_hdf5_multinode( 
         // append a timestamped state to an hdf5 file
         const hid_t outFile_id,
         const int& time_step,
         const double& time,
         const std::vector<double>& phi,
         const std::vector<double>& T,
         const int& Nx_local,
         const std::vector<hsize_t>& dims, // assume dims.size() == 3
         const std::vector<size_t>& idx_start, 
         const std::vector<size_t>& idx_end,
         const hid_t dx_plist_id,
         const int& mynode,
         const int& rootnode,
         const int& totalnodes,
         MPI_Comm comm
         )
{
   int failflag; failflag = 0;
   herr_t status; status = 0;
   hid_t file_space_id, dataset_id, dataspace_id, memspace_id, datasetT_id, dataspaceT_id, memspaceT_id, group_id;
   hid_t attribute_id, attribute_dataspace_id;

   //std::vector<double> output_buffer( Nx_local * dims[1] * dims[2], 0 );
   int ndims; ndims = 3;//dims.size();
   hsize_t offset[3];
   hsize_t offset_local[3];
   hsize_t count[3];
   //hsize_t h5dims[3];
   hsize_t dims_local[3];
   //for ( size_t ii=0; ii < ndims; ++ii) h5dims[ii] = dims[ii];
   dims_local[0] = Nx_local + 2;
   for ( size_t ii=1; ii < ndims; ++ii) dims_local[ii] = dims[ii];

   std::string groupName("Step");
   //std::string groupName("timestep");
   std::string datasetName, datasetNameT;
   std::string attributeName( "Time" );
   hsize_t attribute_dims[1]; attribute_dims[0] = 1;
   
   // convert time_step from into to string
   std::ostringstream sstime_step;
   sstime_step << time_step;
   groupName = "Step" + sstime_step.str();

   datasetName = "/" + groupName + "/phi";
   datasetNameT = "/" + groupName + "/T";

   group_id = H5Gcreate2( outFile_id,
                        groupName.c_str(),
                        H5P_DEFAULT, // hid_t link creation pl_id  
                        H5P_DEFAULT, // hid_t group creation pl_id 
                        H5P_DEFAULT  // hid_t group access pl_id
                           );
   if ( group_id < 0 ) 
   {
      cout << "Error, node " << mynode 
            << " failed to create group " 
            << groupName
            << endl;
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }
   // create time dataspace
   attribute_dataspace_id = H5Screate_simple( 1, attribute_dims, NULL);
   if ( attribute_dataspace_id < 0 ) 
   {
      cout << "Error, node " << mynode 
            << " failed to create attribute dataspace for " 
            << attributeName << " in group "
            << groupName
            << endl;
      H5Sclose( attribute_dataspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // create time attribute
   attribute_id = H5Acreate2(
                           group_id,
                           attributeName.c_str(),
                           H5T_NATIVE_DOUBLE,
                           attribute_dataspace_id,
                           H5P_DEFAULT,
                           H5P_DEFAULT 
                           );
   if ( attribute_id < 0 ) 
   {
      cout << "Error, node " << mynode 
            << " failed to create attribute " 
            << attributeName << " in group "
            << groupName
            << endl;
      H5Aclose( attribute_id );
      H5Sclose( attribute_dataspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // write time attribute to file
   status = H5Awrite( 
                     attribute_id,
                     H5T_NATIVE_DOUBLE,
                     &time
                     );
   if ( status < 0 ) 
   {
      cout << "Error, node " << mynode 
            << " failed to write time attribute to group: "
            << groupName
            << endl;
      H5Aclose( attribute_id );
      H5Sclose( attribute_dataspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   H5Aclose( attribute_id );
   H5Sclose( attribute_dataspace_id );

   //  calculate offset and count
   offset_local[0] = 1; offset_local[1] = 0; offset_local[2] = 0;
   for ( size_t ii=0; ii < ndims; ++ii ) 
   {
      offset[ii] = idx_start[ii];

   }

   count[0] = Nx_local;
   count[1] = dims_local[1];
   count[2] = dims_local[2];

   //hsize_t stride[ndims];
   //for ( size_t ii=0; ii < ndims; ++ii ) stride[ii] = 1;
   //hsize_t block[ndims];
   //for ( size_t ii=0; ii < ndims; ++ii) block[ii] = 1;

   //cout << "node " << mynode << " offset_local: "; // debug
   //for (size_t ii=0; ii < 3 ; ++ii) // debug
   //   cout << offset_local[ii] << " "; // debug
   //cout << endl; // debug
   //cout << "node " << mynode << " count: "; // debug
   //for (size_t ii=0; ii < 3 ; ++ii) // debug
   //   cout << count[ii] << " "; // debug
   //cout << endl; // debug

   // create local memory dataspace
   memspace_id = H5Screate_simple( ndims, dims_local, NULL );
   if ( memspace_id < 0 )
   {
      cout << "Error, node " << mynode 
         << " failed to create memspace_id"  << endl;
      H5Sclose( memspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // create file dataspace
   dataspace_id = H5Screate_simple( ndims, &dims[0], NULL );
   if ( dataspace_id < 0 )
   {
      cout << "Error, node " << mynode 
         << " failed to create dataspace_id" << endl;
      H5Sclose( memspace_id );
      H5Sclose( dataspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }
   
   // create dataset
   dataset_id = H5Dcreate2(
                           outFile_id, 
                           datasetName.c_str(), //"/phi", 
                           H5T_NATIVE_DOUBLE, 
                           dataspace_id, 
                           H5P_DEFAULT, // link creation property list
                           H5P_DEFAULT, // dataset creation property list
                           H5P_DEFAULT // dataset access property list
                           );
   if ( dataset_id < 0 ) 
   {
         cout << "Error, node " << mynode 
            << " failed to create dataset_id from dataspace_id " 
            << endl;
      H5Sclose( memspace_id );
      H5Sclose( dataspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }
   
   // select hyperslab in memory
   status = H5Sselect_hyperslab( 
                                 memspace_id,
                                 H5S_SELECT_SET, // H5S_seloper_t op
                                 offset_local, // const hsize_t *start,
                                 NULL, // identical to having stride == 1
                                 //stride, // const hsize_t *stride,
                                 count, // const hsize_t *count,
                                 NULL // identical to having block==1
                                 //block // const hsize_t *block,
                                 );
   if ( status < 0 ) 
   {
      cout << "Error node " << mynode 
         << " failed to select hyperslab from memspace_id"
         << endl;
      H5Sclose( memspace_id );
      H5Sclose( dataspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // select hyperslab in file
   status = H5Sselect_hyperslab( 
                                 dataspace_id,
                                 H5S_SELECT_SET, // H5S_seloper_t op
                                 offset, // const hsize_t *start,
                                 NULL, // identical to having stride == 1
                                 //stride, // const hsize_t *stride,
                                 count, // const hsize_t *count,
                                 NULL // identical to having block==1
                                 //block // const hsize_t *block,
                                 );
   if ( status < 0 ) 
   {
      cout << "Error node " << mynode 
         << " failed to select hyperslab from dataspace_id"
         << endl;
      H5Sclose( memspace_id );
      H5Sclose( dataspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // write phi memory hyperslab to file hyperslab
   status = H5Dwrite(
         dataset_id, // dataset_id
         H5T_NATIVE_DOUBLE, // mem_type_id
         memspace_id, // H5S_ALL, // mem_space_id
         dataspace_id, // file_space_id
         dx_plist_id, // xfer_plist_id,
         //H5P_DEFAULT, // xfer_plist_id,
         &phi[0] // buf
         //&output_buffer[0] // buf
         );
   if ( status < 0 ) 
   {
      cout << "Error node " << mynode 
         << " failed to write phi with memspace_id dataspace_id"
            << endl;
      H5Sclose( memspace_id );
      H5Sclose( dataspace_id );
      H5Dclose( dataset_id);
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // success in writing phi, so close its H5 data spaces and set
   H5Sclose( memspace_id );
   H5Sclose( dataspace_id );
   H5Dclose( dataset_id );
   
   // create local memory dataspace
   memspaceT_id = H5Screate_simple( ndims, dims_local, NULL );
   if ( memspaceT_id < 0 )
   {
      cout << "Error, node " << mynode 
         << " failed to create memspaceT_id"  << endl;
      H5Sclose( memspaceT_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // create file dataspace
   dataspaceT_id = H5Screate_simple( ndims, &dims[0], NULL );
   if ( dataspaceT_id < 0 )
   {
      cout << "Error, node " << mynode 
         << " failed to create dataspaceT_id" << endl;
      H5Sclose( memspaceT_id );
      H5Sclose( dataspaceT_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // create dataset
   datasetT_id = H5Dcreate2(
                           outFile_id, 
                           datasetNameT.c_str(),
                           H5T_NATIVE_DOUBLE, 
                           dataspaceT_id, 
                           H5P_DEFAULT, // link creation property list
                           H5P_DEFAULT, // dataset creation property list
                           H5P_DEFAULT // dataset access property list
                           );
   if ( datasetT_id < 0 ) 
   {
         cout << "Error, node " << mynode 
            << " failed to create datasetT_id from dataspaceT_id " 
            << endl;
      H5Sclose( memspaceT_id );
      H5Sclose( dataspaceT_id );
      H5Dclose( datasetT_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // select hyperslab in memory
   status = H5Sselect_hyperslab( 
                                 memspaceT_id,
                                 H5S_SELECT_SET, // H5S_seloper_t op
                                 offset_local, // const hsize_t *start,
                                 NULL, // identical to having stride == 1
                                 //stride, // const hsize_t *stride,
                                 count, // const hsize_t *count,
                                 NULL // identical to having block==1
                                 //block // const hsize_t *block,
                                 );
   if ( status < 0 ) 
   {
      cout << "Error node " << mynode 
         << " failed to select hyperslab from memspaceT_id"
         << endl;
      H5Sclose( memspaceT_id );
      H5Sclose( dataspaceT_id );
      H5Dclose( datasetT_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // select hyperslab in file
   status = H5Sselect_hyperslab( 
                                 dataspaceT_id,
                                 H5S_SELECT_SET, // H5S_seloper_t op
                                 offset, // const hsize_t *start,
                                 NULL, // identical to having stride == 1
                                 //stride, // const hsize_t *stride,
                                 count, // const hsize_t *count,
                                 NULL // identical to having block==1
                                 //block // const hsize_t *block,
                                 );
   if ( status < 0 ) 
   {
      cout << "Error node " << mynode 
         << " failed to select hyperslab from dataspaceT_id"
         << endl;
      H5Sclose( memspaceT_id );
      H5Sclose( dataspaceT_id );
      H5Dclose( datasetT_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // write T memory hyperslab to file hyperslab
   status = H5Dwrite(
         datasetT_id, // dataset_id
         H5T_NATIVE_DOUBLE, // mem_type_id
         memspaceT_id, // H5S_ALL, // mem_space_id
         dataspaceT_id, // file_space_id
         dx_plist_id, // xfer_plist_id,
         //H5P_DEFAULT, // xfer_plist_id,
         &T[0] // buf
         //&output_buffer[0] // buf
         );
   if ( status < 0 ) 
   {
      cout << "Error node " << mynode 
         << " failed to write T with memspaceT_id dataspaceT_id" << endl;
      H5Sclose( memspaceT_id );
      H5Sclose( dataspaceT_id );
      H5Dclose( datasetT_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // success in writing T, so close its H5 data spaces and set
   H5Sclose( memspaceT_id );
   H5Sclose( dataspaceT_id );
   H5Dclose( datasetT_id );

   H5Gclose(group_id);
   return EXIT_SUCCESS;
}
#endif
