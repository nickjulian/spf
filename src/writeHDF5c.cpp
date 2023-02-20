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
// File: writeHDF5c.cpp

#ifndef WRITEHDF5C_CPP
#define WRITEHDF5C_CPP

#include "writeHDF5c.hpp"

int SPF_NS::append_phi_to_hdf5_multinode( 
         // append a timestamped state to an hdf5 file
         const hid_t outFile_id,
         const int& time_step,
         const double& time,
         const std::vector<double>& phi,
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
   hid_t file_space_id, dataset_id, dataspace_id, memspace_id, group_id;
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
   std::string datasetName;
   std::string attributeName( "Time" );
   hsize_t attribute_dims[1]; attribute_dims[0] = 1;
   
   // convert time_step from into to string
   std::ostringstream sstime_step;
   sstime_step << time_step;
   groupName = "Step" + sstime_step.str();

   datasetName = "/" + groupName + "/phi";
   
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

   //offset_local[0] = 0; offset_local[1] = 0; offset_local[2] = 0;
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
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }
   H5Sclose( memspace_id );
   H5Sclose( dataspace_id );
   H5Dclose( dataset_id );
   H5Gclose( group_id );
   return EXIT_SUCCESS;
}

int SPF_NS::append_fields_to_hdf5_multinode( 
         // append a timestamped state to an hdf5 file
         const hid_t outFile_id,
         const int& time_step,
         const double& time,
         const std::vector<double>& phi,
         const std::vector<double>& T,
         const std::vector<double>& conc,
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
   hid_t file_space_id, datasetPhi_id, dataspacePhi_id, memspacePhi_id,
         datasetT_id, dataspaceT_id, memspaceT_id,
         datasetConc_id, dataspaceConc_id, memspaceConc_id,
         group_id;
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
   std::string datasetNamePhi, datasetNameT, datasetNameConc;
   std::string attributeName( "Time" );
   hsize_t attribute_dims[1]; attribute_dims[0] = 1;
   
   // convert time_step from into to string
   std::ostringstream sstime_step;
   sstime_step << time_step;
   groupName = "Step" + sstime_step.str();

   datasetNamePhi = "/" + groupName + "/phi";
   datasetNameT = "/" + groupName + "/T";
   datasetNameConc = "/" + groupName + "/conc";
   
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

   // calculate offset and count
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
   memspacePhi_id = H5Screate_simple( ndims, dims_local, NULL );
   if ( memspacePhi_id < 0 )
   {
      cout << "Error, node " << mynode 
         << " failed to create memspacePhi_id"  << endl;
      H5Sclose( memspacePhi_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // create file dataspace
   dataspacePhi_id = H5Screate_simple( ndims, &dims[0], NULL );
   if ( dataspacePhi_id < 0 )
   {
      cout << "Error, node " << mynode 
         << " failed to create dataspacePhi_id" << endl;
      H5Sclose( memspacePhi_id );
      H5Sclose( dataspacePhi_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // create dataset
   datasetPhi_id = H5Dcreate2(
                           outFile_id, 
                           datasetNamePhi.c_str(), //"/phi", 
                           H5T_NATIVE_DOUBLE, 
                           dataspacePhi_id,
                           H5P_DEFAULT, // link creation property list
                           H5P_DEFAULT, // dataset creation property list
                           H5P_DEFAULT // dataset access property list
                           );
   if ( datasetPhi_id < 0 ) 
   {
         cout << "Error, node " << mynode 
            << " failed to create datasetPhi_id from dataspacePhi_id " 
            << endl;
      H5Sclose( memspacePhi_id );
      H5Sclose( dataspacePhi_id );
      H5Dclose( datasetPhi_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // select hyperslab in memory
   status = H5Sselect_hyperslab( 
                                 memspacePhi_id,
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
         << " failed to select hyperslab from memspacePhi_id"
         << endl;
      H5Sclose( memspacePhi_id );
      H5Sclose( dataspacePhi_id );
      H5Dclose( datasetPhi_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // select hyperslab in file
   status = H5Sselect_hyperslab( 
                                 dataspacePhi_id,
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
         << " failed to select hyperslab from dataspacePhi_id"
         << endl;
      H5Sclose( memspacePhi_id );
      H5Sclose( dataspacePhi_id );
      H5Dclose( datasetPhi_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // write phi memory hyperslab to file hyperslab
   status = H5Dwrite(
         datasetPhi_id, // dataset_id
         H5T_NATIVE_DOUBLE, // mem_type_id
         memspacePhi_id, // H5S_ALL, // mem_space_id
         dataspacePhi_id, // file_space_id
         dx_plist_id, // xfer_plist_id,
         //H5P_DEFAULT, // xfer_plist_id,
         &phi[0] // buf
         //&output_buffer[0] // buf
         );
   if ( status < 0 ) 
   {
      cout << "Error node " << mynode 
         << " failed to write phi with memspace_id  dataspacePhi_id"
         << endl;
      H5Sclose( memspacePhi_id );
      H5Sclose( dataspacePhi_id );
      H5Dclose( datasetPhi_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // success in writing phi, so close its H5 data spaces and set
   H5Sclose( memspacePhi_id );
   H5Sclose( dataspacePhi_id );
   H5Dclose( datasetPhi_id );

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

   // create local memory dataspace
   memspaceConc_id = H5Screate_simple( ndims, dims_local, NULL );
   if ( memspaceConc_id < 0 )
   {
      cout << "Error, node " << mynode 
         << " failed to create memspaceConc_id"  << endl;
      H5Sclose( memspaceConc_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }
   // create file dataspace
   dataspaceConc_id = H5Screate_simple( ndims, &dims[0], NULL );
   if ( dataspaceConc_id < 0 )
   {
      cout << "Error, node " << mynode 
         << " failed to create dataspaceConc_id" << endl;
      H5Sclose( memspaceConc_id );
      H5Sclose( dataspaceConc_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // create dataset
   datasetConc_id = H5Dcreate2(
                           outFile_id, 
                           datasetNameConc.c_str(),
                           H5T_NATIVE_DOUBLE, 
                           dataspaceConc_id, 
                           H5P_DEFAULT, // link creation property list
                           H5P_DEFAULT, // dataset creation property list
                           H5P_DEFAULT // dataset access property list
                           );
   if ( datasetConc_id < 0 ) 
   {
         cout << "Error, node " << mynode 
            << " failed to create datasetConc_id from dataspaceConc_id "
            << endl;
      H5Sclose( memspaceConc_id );
      H5Sclose( dataspaceConc_id );
      H5Dclose( datasetConc_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }
   // select hyperslab in memory
   status = H5Sselect_hyperslab( 
                                 memspaceConc_id,
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
         << " failed to select hyperslab from memspaceConc_id"
         << endl;
      H5Sclose( memspaceConc_id );
      H5Sclose( dataspaceConc_id );
      H5Dclose( datasetConc_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }
   // select hyperslab in file
   status = H5Sselect_hyperslab( 
                                 dataspaceConc_id,
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
         << " failed to select hyperslab from dataspaceConc_id"
         << endl;
      H5Sclose( memspaceConc_id );
      H5Sclose( dataspaceConc_id );
      H5Sclose( datasetConc_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }
   // write conc memory hyperslab to file hyperslab
   status = H5Dwrite(
         datasetConc_id, // dataset_id
         H5T_NATIVE_DOUBLE, // mem_type_id
         memspaceConc_id, // H5S_ALL, // mem_space_id
         dataspaceConc_id, // file_space_id
         dx_plist_id, // xfer_plist_id,
         //H5P_DEFAULT, // xfer_plist_id,
         &conc[0] // buf
         //&output_buffer[0] // buf
         );
   if ( status < 0 ) 
   {
      cout << "Error node " << mynode 
         << " failed to write conc with memspace_id dataspaceConc_id"
         << endl;
      H5Sclose( memspaceConc_id );
      H5Sclose( dataspaceConc_id );
      H5Dclose( datasetConc_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // success in writing conc, so close its H5 data spaces and set
   H5Sclose( memspaceConc_id );
   H5Sclose( dataspaceConc_id );
   H5Dclose( datasetConc_id );

   H5Gclose( group_id );
   return EXIT_SUCCESS;
}

int SPF_NS::write_phi_to_hdf5_multinode( 
         // write a single state to the entire '/phi' dataset
         const hid_t outFile_id,
         const std::vector<double>& phi, 
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
   hid_t file_space_id, dataset_id, dataspace_id, memspace_id;

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
   memspace_id = H5Screate_simple( ndims, dims_local, NULL );
   // create dataspace
   if ( memspace_id < 0 )
   {
      cout << "Error, node " << mynode 
         << " failed to create memspace_id"  << endl;
      H5Sclose( memspace_id );
      return EXIT_FAILURE;
   }

   dataspace_id = H5Screate_simple( ndims, &dims[0], NULL );
   if ( dataspace_id < 0 )
   {
      cout << "Error, node " << mynode 
         << " failed to create dataspace_id" << endl;
      H5Sclose( memspace_id );
      H5Sclose( dataspace_id );
      return EXIT_FAILURE;
   }

   dataset_id = H5Dcreate2( 
                           outFile_id, 
                           "/phi", 
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
      return EXIT_FAILURE;
   }

   //offset_local[0] = 0; offset_local[1] = 0; offset_local[2] = 0;
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

   status = H5Sselect_hyperslab( 
                                 memspace_id,
                                 H5S_SELECT_SET, // H5S_seloper_t op
                                 offset_local, // const hsize_t *start,
                                 NULL, // identical to having stride==1
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
      return EXIT_FAILURE;
   }
   status = H5Sselect_hyperslab( 
                                 dataspace_id,
                                 H5S_SELECT_SET, // H5S_seloper_t op
                                 offset, // const hsize_t *start,
                                 NULL, // identical to having stride==1
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
      return EXIT_FAILURE;
   }
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
      return EXIT_FAILURE;
   }
   H5Sclose( memspace_id );
   H5Sclose( dataspace_id );
   H5Dclose( dataset_id );
   return EXIT_SUCCESS;
}

int SPF_NS::append_phi_to_hdf5_singlenode( 
         const hid_t outFile_id,
         const int& time_step,
         const double& time,
         const std::string& field_name,
         const std::vector<hsize_t>& dims, // assume dims.size() == 3
         const std::vector<double>& phi
         )
{
   if ( phi.size() != dims[0] * dims[1] * dims[2] )
   {
      cout << "Error, field.size() != dims[0] * dims[1] * dims[2]"
         << " in function append_phi_to_hdf5_singlenode()"
         << endl;
      return EXIT_FAILURE;
   }

   int failflag; failflag = 0;
   herr_t status; status = 0;
   hid_t file_space_id, dataset_id, dataspace_id, memspace_id, group_id;
   hid_t attribute_id, attribute_dataspace_id;

   //std::vector<double> output_buffer( Nx_local * dims[1] * dims[2], 0 );
   int ndims; ndims = 3;//dims.size();
   hsize_t offset[3];
   hsize_t offset_local[3];
   hsize_t count[3];
   //hsize_t h5dims[3];
   //hsize_t dims_local[3];
   //for ( size_t ii=0; ii < ndims; ++ii) h5dims[ii] = dims[ii];
   //dims_local[0] = Nx_local + 2;
   //for ( size_t ii=1; ii < ndims; ++ii) dims_local[ii] = dims[ii];

   std::string groupName("Step");
   //std::string groupName("timestep");
   std::string datasetName;
   std::string attributeName( "Time" );
   hsize_t attribute_dims[1]; attribute_dims[0] = 1;// attributes: Time

   // convert time_step from into to string
   std::ostringstream sstime_step;
   sstime_step << time_step;
   groupName = "Step" + sstime_step.str();
   datasetName = "/" + groupName + "/" + field_name;
   
   group_id = H5Gcreate2( outFile_id,
                        groupName.c_str(),
                        H5P_DEFAULT, // hid_t link creation pl_id  
                        H5P_DEFAULT, // hid_t group creation pl_id 
                        H5P_DEFAULT  // hid_t group access pl_id
                           );
   if ( group_id < 0 ) 
   {
      cout << "Error, " 
            << " failed to create group " 
            << groupName
            << endl;
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   attribute_dataspace_id = H5Screate_simple( 1, attribute_dims, NULL);
   if ( attribute_dataspace_id < 0 ) 
   {
      cout << "Error, " 
            << " failed to create attribute dataspace for " 
            << attributeName << " in group "
            << groupName
            << endl;
      H5Sclose( attribute_dataspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }
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
      cout << "Error, " 
            << " failed to create attribute " 
            << attributeName << " in group "
            << groupName
            << endl;
      H5Aclose( attribute_id );
      H5Sclose( attribute_dataspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }
   status = H5Awrite( 
                     attribute_id,
                     H5T_NATIVE_DOUBLE,
                     &time
                     );
   if ( status < 0 ) 
   {
      cout << "Error, " 
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

   // create local memory dataspace
   memspace_id = H5Screate_simple( ndims, &dims[0], NULL );
   if ( memspace_id < 0 )
   {
      cout << "Error, " 
         << " failed to create memspace_id"  << endl;
      H5Sclose( memspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   // create file dataspace
   dataspace_id = H5Screate_simple( ndims, &dims[0], NULL );
   if ( dataspace_id < 0 )
   {
      cout << "Error, " 
         << " failed to create dataspace_id" << endl;
      H5Sclose( memspace_id );
      H5Sclose( dataspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

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
         cout << "Error, " 
            << " failed to create dataset_id from dataspace_id " 
            << endl;
      H5Sclose( memspace_id );
      H5Sclose( dataspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }

   count[0] = dims[0];
   count[1] = dims[1];
   count[2] = dims[2];

   //hsize_t stride[ndims];
   //for ( size_t ii=0; ii < ndims; ++ii ) stride[ii] = 1;
   //hsize_t block[ndims];
   //for ( size_t ii=0; ii < ndims; ++ii) block[ii] = 1;

   status = H5Dwrite(
               dataset_id, // dataset_id
               H5T_NATIVE_DOUBLE, // mem_type_id
               memspace_id, // H5S_ALL, // mem_space_id
               dataspace_id, // file_space_id
               //dx_plist_id, // xfer_plist_id,
               H5P_DEFAULT, // xfer_plist_id,
               &phi[0] // buf
               );
   if ( status < 0 ) 
   {
      cout << "Error " 
         << " failed to write field with append_phi_to_hdf5_singlenode()"
         << endl;
      H5Sclose( memspace_id );
      H5Sclose( dataspace_id );
      H5Gclose( group_id );
      return EXIT_FAILURE;
   }
   H5Sclose( memspace_id );
   H5Sclose( dataspace_id );
   H5Dclose( dataset_id );
   H5Gclose( group_id );
   return EXIT_SUCCESS;
}


int SPF_NS::output_vector_to_hdf5(
      const std::vector<double>& xx,
      const string& outFileName
      )
{
   size_t Nx = xx.size();
   hid_t outFile_id, 
         xx_dataset_id, xx_dataspace_id;
   hsize_t dims[1];
   herr_t status;
   if ( !( outFile_id = H5Fcreate( outFileName.c_str(), 
                        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)
         ))
   {
      cerr << "Failed to create file : " << outFileName << endl;
      return EXIT_FAILURE;
   }

   /* create the dataspaces for the datasets */
   dims[0] = Nx;
   xx_dataspace_id = H5Screate_simple(1, dims, NULL);

   /* create the datasets */
   xx_dataset_id = H5Dcreate2(outFile_id, "/t", H5T_NATIVE_DOUBLE, 
      xx_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   /* write variable data to the file. */
   cout << "writing vector<double> data to : " << outFileName << endl;
   status = H5Dwrite( xx_dataset_id, H5T_NATIVE_DOUBLE, 
                           xx_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &xx[0]);
   if ( status < 0 )
   {
      cerr << "Error putting vector<double> data; status " << status 
         << endl;
      return EXIT_FAILURE;
   }

   /* close the datasets */
   status = H5Dclose( xx_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing xx_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   /* close the dataspace */
   status = H5Sclose( xx_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing xx_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   /* close the HDF5 file */
   status = H5Fclose( outFile_id );
   if ( status < 0)
   {
      cerr << "Error closing file; status: " <<  status << endl;
      return EXIT_FAILURE;
   }
   cout << "success writing to file: " << outFileName << endl;

   return EXIT_SUCCESS;
}

int SPF_NS::output_2vectors_to_hdf5(
      const std::vector<double>& tt,
      const std::vector<double>& xx,
      const string& outFileName
      )
{
   size_t Nt = tt.size();
   size_t Nx = xx.size();
   if ( Nt != Nx )
   {
      cout << "error, there are more time points than states, exiting" 
         << endl;
   }

   hid_t outFile_id; 
   hid_t tt_dataset_id, tt_dataspace_id, xx_dataset_id, xx_dataspace_id;

   hsize_t dims[1];
   //hsize_t dims[2];
   herr_t status;

   if ( !( outFile_id = H5Fcreate( outFileName.c_str(), 
                        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)
         ))
   {
      cerr << "Failed to create file : " << outFileName << endl;
      return EXIT_FAILURE;
   }

   /* create the dataspaces for the datasets */
   dims[0] = Nx;
   //dims[1] = Nt;

   // H5Screate_simple( int rank, const hsize_t* current_dims, 
   //                   const hsize_t* maximum_dims)
   xx_dataspace_id = H5Screate_simple(1, dims, NULL);
   tt_dataspace_id = H5Screate_simple(1, dims, NULL);

   /* create the datasets */
   tt_dataset_id = H5Dcreate2(outFile_id, "/t", H5T_NATIVE_DOUBLE, 
      tt_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   xx_dataset_id = H5Dcreate2(outFile_id, "/x", H5T_NATIVE_DOUBLE, 
      xx_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   /* write variable data to the file. */
   cout << "writing time data to : " << outFileName << endl;
   status = H5Dwrite( tt_dataset_id, H5T_NATIVE_DOUBLE, 
                           tt_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &tt[0]);
   if ( status < 0 )
   {
      cerr << "Error putting time data; status " << status 
         << endl;
      return EXIT_FAILURE;
   }

   cout << "writing state data to : " << outFileName << endl;
   status = H5Dwrite( xx_dataset_id, H5T_NATIVE_DOUBLE, 
                           xx_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &xx[0]);
   if ( status < 0 )
   {
      cerr << "Error putting time data; status " << status 
         << endl;
      return EXIT_FAILURE;
   }

   /* close the datasets */
   status = H5Dclose( tt_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing tt_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dclose( xx_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing xx_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   /* close the dataspace */
   status = H5Sclose( tt_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing tt_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Sclose( xx_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing xx_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   /* close the HDF5 file */
   status = H5Fclose( outFile_id );
   if ( status < 0)
   {
      cerr << "Error closing file; status: " <<  status << endl;
      return EXIT_FAILURE;
   }
   cout << "success writing to file: " << outFileName << endl;

   return EXIT_SUCCESS;
}

int SPF_NS::output_path1D_to_hdf5(
      const std::vector<double>& xx,
      const std::vector<double>& tt,
      const string& outFilePrefix
      )
{
   /* prepare input data for writing */

   size_t Nx = xx.size();
   size_t Nt = tt.size();
   size_t Npoints;

   if ( (Nx == Nt) )
      Npoints = Nx;
   else
   {
      cerr << "Error: Nx != Nt" << endl;
      return  EXIT_FAILURE;
   }

   string outFileName = outFilePrefix + "_1D.h5";

   /* open hdf5 file */
   hid_t outFile_id, 
         xx_dataset_id, xx_dataspace_id,
         tt_dataset_id, tt_dataspace_id;
   hsize_t dims[1];
   herr_t status;

   if ( !( outFile_id = H5Fcreate( outFileName.c_str(), 
                        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)
         ))
   {
      cerr << "Failed to create file : " << outFileName << endl;
      return EXIT_FAILURE;
   }

   /* create the dataspaces for the datasets */
   dims[0] = Npoints;

   xx_dataspace_id = H5Screate_simple(1, dims, NULL);
   tt_dataspace_id = H5Screate_simple(1, dims, NULL);

   /* create the datasets */
   xx_dataset_id = H5Dcreate2(outFile_id, "/xx1D", H5T_NATIVE_DOUBLE, 
      xx_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   tt_dataset_id = H5Dcreate2(outFile_id, "/tt1D", H5T_NATIVE_DOUBLE, 
      tt_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   /* write variable data to the file. */
   cout << "writing 1-D path data to : " << outFileName << endl;
   status = H5Dwrite( xx_dataset_id, H5T_NATIVE_DOUBLE, 
                           xx_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &xx[0]);
   if ( status < 0 )
   {
      cerr << "Error putting 1-D path xx data; status " << status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dwrite( tt_dataset_id, H5T_NATIVE_DOUBLE, 
                           tt_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &tt[0]);
   if ( status < 0 )
   {
      cerr << "Error putting 1-D path tt data; status " << status << endl;
      return EXIT_FAILURE;
   }

   /* close the datasets */
   status = H5Dclose( xx_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing xx_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dclose( tt_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing tt_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }

   /* close the dataspace */
   status = H5Sclose( xx_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing xx_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Sclose( tt_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing tt_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }

   /* close the HDF5 file */
   status = H5Fclose( outFile_id );
   if ( status < 0)
   {
      cerr << "Error closing file; status: " <<  status << endl;
      return EXIT_FAILURE;
   }

   cout << "success writing to file: " << outFileName << endl;

   return EXIT_SUCCESS;
}

int SPF_NS::output_path2D_to_hdf5(
      const std::vector<double>& xx,
      const std::vector<double>& yy,
      const std::vector<double>& tt,
      const string& outFilePrefix
      )
{
   /* prepare input data for writing */

   size_t Nx = xx.size();
   size_t Ny = yy.size();
   size_t Nt = tt.size();
   size_t Npoints;

   if ( (Nx == Nt) && (Ny == Nt) )
      Npoints = Nx;
   else
   {
      cerr << "Error: Nx != Nt or Ny != Nt" << endl;
      return  EXIT_FAILURE;
   }

   string outFileName = outFilePrefix + "_2D.h5";

   /* open hdf5 file */
   hid_t outFile_id, 
         xx_dataset_id, xx_dataspace_id,
         yy_dataset_id, yy_dataspace_id,
         tt_dataset_id, tt_dataspace_id;
   hsize_t dims[1];
   herr_t status;

   if ( !( outFile_id = H5Fcreate( outFileName.c_str(), 
                        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)
         ))
   {
      cerr << "Failed to create file : " << outFileName << endl;
      return EXIT_FAILURE;
   }

   /* create the dataspaces for the datasets */
   dims[0] = Npoints;

   xx_dataspace_id = H5Screate_simple(1, dims, NULL);
   yy_dataspace_id = H5Screate_simple(1, dims, NULL);
   tt_dataspace_id = H5Screate_simple(1, dims, NULL);

   /* create the datasets */
   xx_dataset_id = H5Dcreate2(outFile_id, "/xx2D", H5T_NATIVE_DOUBLE, 
      xx_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   yy_dataset_id = H5Dcreate2(outFile_id, "/yy2D", H5T_NATIVE_DOUBLE, 
      yy_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   tt_dataset_id = H5Dcreate2(outFile_id, "/tt2D", H5T_NATIVE_DOUBLE, 
      tt_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   ///* write variable data to the file. */
   cout << "writing 2-D path data to : " << outFileName << endl;
   status = H5Dwrite( xx_dataset_id, H5T_NATIVE_DOUBLE, 
                           xx_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &xx[0]);
   if ( status < 0 )
   {
      cerr << "Error putting 2-D path xx data; status " << status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dwrite( yy_dataset_id, H5T_NATIVE_DOUBLE, 
                           yy_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &yy[0]);
   if ( status < 0 )
   {
      cerr << "Error putting 2-D path yy data; status " << status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dwrite( tt_dataset_id, H5T_NATIVE_DOUBLE, 
                           tt_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &tt[0]);
   if ( status < 0 )
   {
      cerr << "Error putting 2-D path tt data; status " << status << endl;
      return EXIT_FAILURE;
   }

   /* close the datasets */
   status = H5Dclose( xx_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing xx_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dclose( yy_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing yy_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dclose( tt_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing tt_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }

   /* close the dataspace */
   status = H5Sclose( xx_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing xx_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Sclose( yy_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing yy_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Sclose( tt_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing tt_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }

   /* close the HDF5 file */
   status = H5Fclose( outFile_id );
   if ( status < 0)
   {
      cerr << "Error closing file; status: " <<  status << endl;
      return EXIT_FAILURE;
   }

   cout << "success writing to file: " << outFileName << endl;

   return EXIT_SUCCESS;
}

int SPF_NS::output_path3D_to_hdf5(
      const std::vector<double>& xx,
      const std::vector<double>& yy,
      const std::vector<double>& zz,
      const std::vector<double>& tt,
      const string& outFilePrefix
      )
{
   /* prepare input data for writing */

   size_t Nx = xx.size();
   size_t Ny = yy.size();
   size_t Nz = zz.size();
   size_t Nt = tt.size();
   size_t Npoints;

   if ( (Nx == Nt) && (Ny == Nt) && (Nz == Nt))
      Npoints = Nx;
   else
   {
      cerr << "Error: Nx != Nt or Ny != Nt or Nz != Nt" << endl;
      return  EXIT_FAILURE;
   }

   string outFileName = outFilePrefix + "_3D.h5";

   /* open hdf5 file */
   hid_t outFile_id, 
         xx_dataset_id, xx_dataspace_id,
         yy_dataset_id, yy_dataspace_id,
         zz_dataset_id, zz_dataspace_id,
         tt_dataset_id, tt_dataspace_id;
   hsize_t dims[1];
   herr_t status;

   if ( !( outFile_id = H5Fcreate( outFileName.c_str(), 
                        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)
         ))
   {
      cerr << "Failed to create file : " << outFileName << endl;
      return EXIT_FAILURE;
   }

   /* create the dataspaces for the datasets */
   dims[0] = Npoints;

   xx_dataspace_id = H5Screate_simple(1, dims, NULL);
   yy_dataspace_id = H5Screate_simple(1, dims, NULL);
   zz_dataspace_id = H5Screate_simple(1, dims, NULL);
   tt_dataspace_id = H5Screate_simple(1, dims, NULL);

   /* create the datasets */
   xx_dataset_id = H5Dcreate2(outFile_id, "/xx2D", H5T_NATIVE_DOUBLE, 
      xx_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   yy_dataset_id = H5Dcreate2(outFile_id, "/yy2D", H5T_NATIVE_DOUBLE, 
      yy_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   zz_dataset_id = H5Dcreate2(outFile_id, "/zz2D", H5T_NATIVE_DOUBLE, 
      zz_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   tt_dataset_id = H5Dcreate2(outFile_id, "/tt2D", H5T_NATIVE_DOUBLE, 
      tt_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   ///* write variable data to the file. */
   cout << "writing 3-D path data to : " << outFileName << endl;
   status = H5Dwrite( xx_dataset_id, H5T_NATIVE_DOUBLE, 
                           xx_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &xx[0]);
   if ( status < 0 )
   {
      cerr << "Error putting 2-D path xx data; status " << status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dwrite( yy_dataset_id, H5T_NATIVE_DOUBLE, 
                           yy_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &yy[0]);
   if ( status < 0 )
   {
      cerr << "Error putting 2-D path yy data; status " << status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dwrite( zz_dataset_id, H5T_NATIVE_DOUBLE, 
                           zz_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &zz[0]);
   if ( status < 0 )
   {
      cerr << "Error putting 2-D path zz data; status " << status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dwrite( tt_dataset_id, H5T_NATIVE_DOUBLE, 
                           tt_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &tt[0]);
   if ( status < 0 )
   {
      cerr << "Error putting 2-D path tt data; status " << status << endl;
      return EXIT_FAILURE;
   }

   /* close the datasets */
   status = H5Dclose( xx_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing xx_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dclose( yy_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing yy_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dclose( zz_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing zz_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dclose( tt_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing tt_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }

   /* close the dataspace */
   status = H5Sclose( xx_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing xx_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Sclose( yy_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing yy_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Sclose( zz_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing zz_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Sclose( tt_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing tt_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }

   /* close the HDF5 file */
   status = H5Fclose( outFile_id );
   if ( status < 0)
   {
      cerr << "Error closing file; status: " <<  status << endl;
      return EXIT_FAILURE;
   }

   cout << "success writing to file: " << outFileName << endl;

   return EXIT_SUCCESS;
}

int SPF_NS::output_odf_to_hdf5(
      const double* const odfRaw,
      const std::vector<double>& omegaIn, // omega domain
      const std::vector<double>& thetaIn, // theta domain
      const std::vector<double>& phiIn, // phi domain
      const string& outFilePrefix
      )
{
   /* prepare input data for writing */

   size_t Nomega = omegaIn.size();
   size_t Ntheta = thetaIn.size();
   size_t Nphi = phiIn.size();
   size_t Npoints;

   if ( (Nomega == Ntheta) && (Nomega == Nphi))
      Npoints = Nomega;
   else
   {
      cerr << "Error: Nomega != Ntheta or Nomega != Nphi" << endl;
      return  EXIT_FAILURE;
   }

   //double omegaDomain[ Npoints ];
   //double thetaDomain[ Npoints ];
   //double phiDomain[ Npoints ];

   //std::vector<double>::const_iterator omegaItr = omegaIn.begin();
   //std::vector<double>::const_iterator thetaItr = thetaIn.begin();
   //std::vector<double>::const_iterator phiItr = phiIn.begin();

   //for ( size_t i=0; i< Nomega; ++i)
   //{
   //   omegaDomain[i] = *omegaItr;
   //   ++omegaItr;
   //}
   //for ( size_t i=0; i< Ntheta; ++i)
   //{
   //   thetaDomain[i] = *thetaItr;
   //   ++thetaItr;
   //}
   //for ( size_t i=0; i< Nphi; ++i)
   //{
   //   phiDomain[i] = *phiItr;
   //   ++phiItr;
   //}
   //for ( size_t i=0; i< Npoints; ++i)
   //{
   //   omegaDomain[i] = *omegaItr;
   //   ++omegaItr;
   //   thetaDomain[i] = *thetaItr;
   //   ++thetaItr;
   //   phiDomain[i] = *phiItr;
   //   ++phiItr;
   //}

   string outFileName = outFilePrefix + "_odf.h5";

   /* open hdf5 file */
   hid_t outFile_id, odf_dataset_id, odf_dataspace_id,
         omega_dataset_id, omega_dataspace_id,
         theta_dataset_id, theta_dataspace_id,
         phi_dataset_id, phi_dataspace_id;
   hsize_t dims[1];
   //hsize_t dims[2];
   herr_t status;

   //int odf_varid;
   //int omega_varid, theta_varid, phi_varid;
   //int status; // returned status of H5F functions

   if ( !(outFile_id = H5Fcreate( outFileName.c_str(), 
                        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)
         //status = nc_create( outFileName.c_str(), 
         //   NC_NETCDF4|NC_CLOBBER, &ncid) 
         ))
   {
      cerr << "Failed to create file : " << outFileName << endl;
      //cerr << nc_strerror( status ) << endl;
      //cerr <<  status << endl;
      return EXIT_FAILURE;
   }
   /* create the dataspaces for the datasets */
   dims[0] = Npoints;
   //dims[1] = 1;//Npoints;
   //dims[2] = Npoints;
   //dims[3] = Npoints;
   //if (!( odf_dataspace_id = H5Screate_simple(1, dims, NULL)))
   //      {
   //   cerr << "Failed to create odf_dataspace; id " << odf_dataspace_id << endl;
   //      }
   odf_dataspace_id = H5Screate_simple(1, dims, NULL);
   omega_dataspace_id = H5Screate_simple(1, dims, NULL);
   theta_dataspace_id = H5Screate_simple(1, dims, NULL);
   phi_dataspace_id = H5Screate_simple(1, dims, NULL);

   /* create the datasets */
   odf_dataset_id = H5Dcreate2(outFile_id, "/odf", H5T_NATIVE_DOUBLE, 
      odf_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   omega_dataset_id = H5Dcreate2(outFile_id, "/omega", H5T_NATIVE_DOUBLE, 
      omega_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   theta_dataset_id = H5Dcreate2(outFile_id, "/theta", H5T_NATIVE_DOUBLE, 
      theta_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   phi_dataset_id = H5Dcreate2(outFile_id, "/phi", H5T_NATIVE_DOUBLE, 
      phi_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   /* write variable data to the file. */
   cout << "writing ODF data to : " << outFileName << endl;
   status = H5Dwrite( odf_dataset_id, H5T_NATIVE_DOUBLE, 
                           //H5S_ALL, outFile_id, H5P_DEFAULT, 
                           odf_dataspace_id, 
                           //outFile_id, 
                           H5S_ALL,
                           H5P_DEFAULT, 
                           &odfRaw[0]);
   if ( status < 0 )
   {
      cerr << "Error putting odf; status " << status << endl;
      return EXIT_FAILURE;
   }
   cout << "writing omega data to : " << outFileName << endl;
   status = H5Dwrite( omega_dataset_id, H5T_NATIVE_DOUBLE, 
                           omega_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &omegaIn[0]);
   if ( status < 0 )
   {
      cerr << "Error putting omega; status " << status << endl;
      return EXIT_FAILURE;
   }
   cout << "writing theta data to : " << outFileName << endl;
   status = H5Dwrite( theta_dataset_id, H5T_NATIVE_DOUBLE, 
                           theta_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &thetaIn[0]);
   if ( status < 0 )
   {
      cerr << "Error putting theta; status " << status << endl;
      return EXIT_FAILURE;
   }
   cout << "writing phi data to : " << outFileName << endl;
   status = H5Dwrite( phi_dataset_id, H5T_NATIVE_DOUBLE, 
                           phi_dataspace_id, H5S_ALL, H5P_DEFAULT, 
                           &phiIn[0]);
   if ( status < 0 )
   {
      cerr << "Error putting phi; status " << status << endl;
      return EXIT_FAILURE;
   }

   /* close the datasets */
   status = H5Dclose( odf_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing odf_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dclose( omega_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing omega_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dclose( theta_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing theta_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Dclose( phi_dataset_id );
   if ( status < 0)
   {
      cerr << "Error closing phi_dataset; status " <<  status << endl;
      return EXIT_FAILURE;
   }

   /* close the dataspace */
   status = H5Sclose( odf_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing odf_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Sclose( omega_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing omega_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Sclose( theta_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing theta_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }
   status = H5Sclose( phi_dataspace_id );
   if ( status < 0)
   {
      cerr << "Error closing phi_dataspace; status " <<  status << endl;
      return EXIT_FAILURE;
   }

   /* close the HDF5 file */
   status = H5Fclose( outFile_id );
   if ( status < 0)
   {
      cerr << "Error closing file; status: " <<  status << endl;
      return EXIT_FAILURE;
   }

   cout << "success writing to file: " << outFileName << endl;

   return EXIT_SUCCESS;
}

#endif
