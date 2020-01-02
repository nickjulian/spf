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
      int& write_period,
      int_flags& flags,
      string& output_prefix,
      string& input_field_name
      )
{
   return read_cmdline_options(
         args,
         dt,
         Nt,
         write_period,
         flags,
         output_prefix,
         input_field_name,
         0, 0,
         MPI_COMM_WORLD);
}

int SPF_NS::read_cmdline_options(
      const std::vector<string>& args,
      double& dt,
      int& Nt,
      int& write_period,
      int_flags& flags,
      string& output_prefix,
      string& input_field_name,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   int Nv; Nv = 1;
   return read_cmdline_options(
         args,
         dt,
         Nt,
         Nv,
         write_period,
         flags,
         output_prefix,
         input_field_name,
         0, 0,
         MPI_COMM_WORLD);
}

int SPF_NS::read_cmdline_options(
      const std::vector<string>& args,
      double& dt,
      int& Nt,
      int& Nv,
      int& write_period,
      int_flags& flags,
      string& output_prefix,
      string& input_field_name,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   string parameter_filename;

   // order of reading parameters allows cmdline to override those in file
   for ( size_t idx=1; idx < args.size(); idx++)
   {
      //std::cout << "reading argument : " << idx  // debug
      //   << " " << args[idx] << std::endl; // debug
      if ( args[idx] == "--parameter_file" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> parameter_filename;
            if ( read_parameter_file(
                     parameter_filename,
                     flags,
                     dt,
                     Nt,
                     Nv,
                     write_period,
                     output_prefix,
                     input_field_name,
                     mynode,
                     rootnode,
                     comm
                     ) != EXIT_SUCCESS)
            {
               cout << "Failed to read parameter file : "
                  << parameter_filename << " . Exiting" << endl;
               return EXIT_FAILURE;
            }
            flags.parameter_file = 1;
            idx += 1;
         }
         //else
         //   if ( mynode == rootnode )
         //   {
         //      std::cout 
         //         << "missing argument to option '--parameter_file'" 
         //         << std::endl;
         //   }
      }
      else if ( args[idx] == "-o" )
      {
         if ( idx + 1 < args.size()) 
         {
            output_prefix = args[idx + 1];
            flags.output_prefix = 1;   // true
            idx += 1;
         }
         else
            if ( mynode == rootnode )
            {
               std::cout << "missing argument to option '-o'" << std::endl;
            }
      }
      else if ( args[idx] == "-i" )
      {
         if ( idx + 1 < args.size()) 
         {
            input_field_name = string( args[idx + 1] );
            flags.input_field = 1;  // true
            idx += 1;
         }
      }
      else if ( args[idx] == "-dt" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> dt;
            flags.dt = 1;
            idx += 1;
         }
      }
      else if ( args[idx] == "-Nt" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> Nt;
            flags.Nt = 1;
            idx += 1;
         }
      }
      else if ( args[idx] == "-Nv" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> Nv;
            flags.Nv = 1;
            idx += 1;
         }
      }
      else if ( args[idx] == "-wp" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> write_period;
            flags.wp = 1;
            idx += 1;
         }
      }
      else if ( args[idx] == "-stat" ) 
      {
         flags.calcstat= 1;
      }
      else if ( args[idx] == "--debug" ) 
      {
         flags.debug = 1;
      }
   }

   if ( !((flags.output_prefix == 1) & (flags.input_field == 1) ))
   {
      if ( mynode == rootnode )
         cout << "Error: insufficient arguments;" 
            << " require both input field (-i <file>)"
            << " and output prefix (-o <path/file> )" << endl;

      return EXIT_FAILURE;
   }
   if ( flags.dt == 0 )
   {
      if ( mynode == rootnode )
         cout << "warning: time increment (-dt <#>) not specified,"
           << " assigning dt = 1"
            << endl;
      dt = 1.0;
   }
   if ( flags.Nt == 0 )
   {
      if ( mynode == rootnode )
         cout << "Error: number of time steps (-Nt <#> ) not specified"
            << endl;
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
} // read_cmdline_options()


#endif
