/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: read_cmdline.cpp

#ifndef READ_CMDLINE_CPP
#define READ_CMDLINE_CPP

#include "read_cmdline.hpp"

int SPF_NS::read_cmdline_options(
      const std::vector<string>& args,
      double& dt,
      int& Nt,
      double& rate_scale_factor,
      int& write_period,
      bool& flag_calcstat,
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
   bool flag_dt; flag_dt = false;
   bool flag_Nt; flag_Nt = false;
   bool flag_r; flag_r = false;
   bool flag_wp; flag_wp = false;

   for ( size_t idx=1; idx < args.size(); idx++)
   {
      if ( args[idx] == "-o" )
         if ( idx + 1 < args.size()) 
         {
            output_prefix = args[idx + 1];
            flag_outprefix = true;
            idx += 1;
         }
      if ( args[idx] == "-i" )
         if ( idx + 1 < args.size()) 
         {
            input_field_name = string( args[idx + 1] );
            flag_input_field = true;
            idx += 1;
         }
      if ( args[idx] == "-dt" )
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> dt;
            flag_dt = true;
         }
      if ( args[idx] == "-Nt" )
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> Nt;
            flag_Nt = true;
         }
      if ( args[idx] == "-r" )
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> rate_scale_factor;
            flag_r = true;
         }
      if ( args[idx] == "-wp" )
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> write_period;
            flag_wp = true;
         }
      if ( args[idx] == "-stat" ) flag_calcstat= true;
   }
   if ( !(flag_outprefix & flag_input_field ))
   {
      if ( mynode == rootnode )
         cout << "Error: insufficient arguments ..." << endl;

      return EXIT_FAILURE;
   }
   if ( !flag_dt )
   {
      if ( mynode == rootnode )
         cout << "warning: time increment (-dt <#>) not specified,"
           << " assigning dt = 1"
            << endl;
      dt = 1.0;
   }
   if ( !flag_Nt )
   {
      if ( mynode == rootnode )
         cout << "Error: number of time steps (-Nt <#> ) not specified"
            << endl;
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
} // read_cmdline_options()

#endif
