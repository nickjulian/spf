/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: read_cmdline.cpp

#ifndef READ_CMDLINE_CPP
#define READ_CMDLINE_CPP

#include "read_cmdline.hpp"

int SPF_NS::read_cmdline_options(
      const std::vector<string>& args,
      string& output_prefix,
      string& input_field_name,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   // define flags to check input requirements
   bool flag_outprefix; flag_outprefix= false;
   bool flag_input_field; flag_input_field = false;

   for ( size_t idx=1; idx < args.size(); idx++)
   {
      if ( args[idx] == "-o" )
         if ( idx + 1 < args.size()) 
         {
            output_prefix = args[idx + 1];
            flag_outprefix = true;
         }
      if ( args[idx] == "-i" )
         if ( idx + 1 < args.size()) 
         {
            input_field_name = args[idx + 1];
            flag_input_field = true;
         }
   }
   if ( !(flag_outprefix & flag_input_field ))
   {
      if ( mynode == rootnode )
         cout << "Error: insufficient arguments ..." << endl;

      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
} // read_cmdline_options()

#endif
