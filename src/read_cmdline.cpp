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
   string datasetPath = "/phi";
   return read_cmdline_options(
         args,
         dt,
         Nt,
         write_period,
         flags,
         output_prefix,
         input_field_name,
         datasetPath,
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
      string& datasetPath
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
         datasetPath,
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
   string datasetPath = "/phi";
   return read_cmdline_options(
         args,
         dt,
         Nt,
         Nv,
         write_period,
         flags,
         output_prefix,
         input_field_name,
         datasetPath,
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
      string& datasetPath,
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
         datasetPath,
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
      string& datasetPath,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   // variables created just to satisfy the full function prototype
   double hh_x;
   double ww;
   double shape_constant;
   double mobility;
   double kappa;
   double c_alpha;
   double c_beta;

   return read_cmdline_options(
         args,
         dt,
         Nt,
         Nv,
         hh_x,
         ww,
         shape_constant,
         mobility,
         kappa,
         c_alpha,
         c_beta,
         write_period,
         flags,
         output_prefix,
         input_field_name,
         datasetPath,
         mynode, rootnode,
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
   // variables created just to satisfy the full function prototype
   double hh_x;
   double ww;
   double shape_constant;
   double mobility;
   double kappa;
   double c_alpha;
   double c_beta;
   string datasetPath = "/phi";

   return read_cmdline_options(
         args,
         dt,
         Nt,
         Nv,
         hh_x,
         ww,
         shape_constant,
         mobility,
         kappa,
         c_alpha,
         c_beta,
         write_period,
         flags,
         output_prefix,
         input_field_name,
         datasetPath,
         mynode, rootnode,
         MPI_COMM_WORLD);
}

int SPF_NS::read_cmdline_options(
      const std::vector<string>& args,
      double& dt,
      int& Nt,
      int& Nv,
      double& hh_x,
      double& ww,
      double& shape_constant,
      double& mobility,
      double& kappa,
      double& c_alpha,
      double& c_beta,
      int& write_period,
      int_flags& flags,
      string& output_prefix,
      string& input_field_name,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   string datasetPath = "/phi";

   return read_cmdline_options(
         args,
         dt,
         Nt,
         Nv,
         hh_x,
         ww,
         shape_constant,
         mobility,
         kappa,
         c_alpha,
         c_beta,
         write_period,
         flags,
         output_prefix,
         input_field_name,
         datasetPath,
         mynode, rootnode,
         MPI_COMM_WORLD);
}

int SPF_NS::read_cmdline_options(
      const std::vector<string>& args,
      double& dt,
      int& Nt,
      int& Nv,
      double& hh_x,
      double& ww,
      double& shape_constant,
      double& mobility,
      double& kappa,
      double& c_alpha,
      double& c_beta,
      int& write_period,
      int_flags& flags,
      string& output_prefix,
      string& input_field_name,
      string& datasetPath,
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
                     hh_x,
                     ww,
                     shape_constant,
                     mobility,
                     kappa,
                     c_alpha,
                     c_beta,
                     write_period,
                     output_prefix,
                     input_field_name,
                     datasetPath,
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
      else if ( args[idx] == "-datasetPath" )
      {
         if ( idx + 1 < args.size()) 
         {
            datasetPath = string( args[idx + 1] );
            flags.datasetPath = 1;  // true
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
      else if ( args[idx] == "-debug" ) 
      {
         flags.debug = 1;
      }
      else if ( args[idx] == "-mesh-size" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> hh_x;
            idx += 1;
         }
      }
      else if ( args[idx] == "-order-energy" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> ww;
            idx += 1;
         }
      }
      else if ( args[idx] == "-shape-constant" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> shape_constant;
            idx += 1;
         }
      }
      else if ( args[idx] == "-mobility" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> mobility;
            idx += 1;
         }
      }
      else if ( args[idx] == "-kappa" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> kappa;
            idx += 1;
         }
      }
      else if ( args[idx] == "-c-alpha" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> c_alpha;
            idx += 1;
         }
      }
      else if ( args[idx] == "-c-beta" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> c_beta;
            idx += 1;
         }
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

int SPF_NS::read_cmdline_options(
      const std::vector<string>& args,
      double& dt,
      int& Nt,
      int& Nv,
      double& hh_x,
      double& c1,
      double& c2,
      double& c3,
      double& c4,
      double& c5,
      double& c6,
      double& cPrefactor,
      double& alpha,
      double& T0,
      double& cbase,
      double& LL,
      double& orderEnergy,
      double& D_T,
      double& M_phi,
      double& M_conc,
      std::vector<double>& fieldValueLimits,
      int& write_period,
      int_flags& flags,
      string& output_prefix,
      string& input_field_file_name,
      string& datasetPathPhi,
      string& datasetPathT,
      string& datasetPathConc,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   string parameter_filename;

   if ( fieldValueLimits.size() != 5)
   {
      std::cout << "Error: "
         "SPF_NS::read_cmdline_options() fieldValueLimits.size() != 5"
         << std::endl;
      return EXIT_FAILURE;
   }

   // order of reading parameters allows cmdline to override those in file
   for ( size_t idx=1; idx < args.size(); idx++)
   {
      //std::cout << "reading argument : " << idx  // debug
      //   << " " << args[idx] << std::endl; // debug
      std::string datasetPath;
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
                     hh_x,
                     c1,
                     c2,
                     c3,
                     c4,
                     c5,
                     c6,
                     cPrefactor,
                     alpha,
                     T0,
                     cbase,
                     LL,
                     orderEnergy,
                     D_T,
                     M_phi,
                     M_conc,
                     fieldValueLimits,
                     write_period,
                     output_prefix,
                     input_field_file_name,
                     datasetPath,
                     datasetPathPhi,
                     datasetPathT,
                     datasetPathConc,
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
               std::cout << "missing argument to option '-o'"
                  << std::endl;
            }
      }
      else if ( args[idx] == "-i" )
      {
         if ( idx + 1 < args.size()) 
         {
            input_field_file_name = string( args[idx + 1] );
            flags.input_field = 1;  // true
            idx += 1;
         }
      }
      else if ( args[idx] == "-datasetPathPhi" )
      {
         if ( idx + 1 < args.size()) 
         {
            datasetPathPhi = string( args[idx + 1] );
            flags.datasetPathPhi = 1;  // true
            idx += 1;
         }
      }
      else if ( args[idx] == "-datasetPathT" )
      {
         if ( idx + 1 < args.size()) 
         {
            datasetPathT = string( args[idx + 1] );
            flags.datasetPathT = 1;  // true
            idx += 1;
         }
      }
      else if ( args[idx] == "-datasetPathConc" )
      {
         if ( idx + 1 < args.size()) 
         {
            datasetPathConc= string( args[idx + 1] );
            flags.datasetPathConc= 1;  // true
            idx += 1;
         }
      }
      else if ( args[idx] == "-datasetPath" )
      {
         if ( idx + 1 < args.size()) 
         {
            datasetPath = string( args[idx + 1] );
            flags.datasetPath = 1;  // true
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
      else if ( args[idx] == "-write_period" )
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
      else if ( args[idx] == "-debug" ) 
      {
         flags.debug = 1;
      }
      else if ( args[idx] == "-mesh-size" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> hh_x;
            idx += 1;
         }
      }
      else if ( args[idx] == "-orderEnergy" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> orderEnergy;
            idx += 1;
         }
      }
      else if ( args[idx] == "-c1" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> c1;
            idx += 1;
         }
      }
      else if ( args[idx] == "-c2" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> c2;
            idx += 1;
         }
      }
      else if ( args[idx] == "-c3" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> c3;
            idx += 1;
         }
      }
      else if ( args[idx] == "-c4" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> c4;
            idx += 1;
         }
      }
      else if ( args[idx] == "-c5" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> c5;
            idx += 1;
         }
      }
      else if ( args[idx] == "-c6" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> c6;
            idx += 1;
         }
      }
      else if ( args[idx] == "-cPrefactor" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> cPrefactor;
            idx += 1;
         }
      }
      else if ( args[idx] == "-alpha" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> alpha;
            idx += 1;
         }
      }
      else if ( args[idx] == "-T0" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> T0;
            idx += 1;
         }
      }
      else if ( args[idx] == "-cbase" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> cbase;
            idx += 1;
         }
      }
      else if ( args[idx] == "-LL" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> LL;
            idx += 1;
         }
      }
      else if ( args[idx] == "-D_T" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> D_T;
            idx += 1;
         }
      }
      else if ( args[idx] == "-M_phi" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> M_phi;
            idx += 1;
         }
      }
      else if ( args[idx] == "-M_conc" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> M_conc;
            idx += 1;
         }
      }
      else if ( args[idx] == "-phi_lower_limit" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> fieldValueLimits[0];
            idx += 1;
         }
      }
      else if ( args[idx] == "-phi_upper_limit" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> fieldValueLimits[1];
            idx += 1;
         }
      }
      else if ( args[idx] == "-conc_lower_limit" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> fieldValueLimits[2];
            idx += 1;
         }
      }
      else if ( args[idx] == "-conc_upper_limit" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> fieldValueLimits[3];
            idx += 1;
         }
      }
      else if ( args[idx] == "-T_lower_limit" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> fieldValueLimits[4];
            idx += 1;
         }
      }
      //else if ( args[idx] == "-T_upper_limit" )
      //{
      //   if ( idx + 1 < args.size()) 
      //   {
      //      istringstream( args[idx + 1] ) >> fieldValueLimits[5];
      //      idx += 1;
      //   }
      //}
   }

   if ((flags.datasetPathPhi == 0) || (flags.datasetPathT == 0)
            || (flags.datasetPathConc == 0))

   {
      if ( mynode == rootnode )
         cout << "Error: insufficient arguments;" 
            << " require dataset paths for three fields Phi, T, Conc"
            << " (-datasetPathPhi </phi>,"
           << " -datasetPathT </T>, -datasetPathConc </conc>)"
            << " and output prefix (-o <path/file> )" << endl;
      return EXIT_FAILURE;
   }
   else 
   {
      flags.datasetPath = 1; // in case it's needed by subsequent functions
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

int SPF_NS::read_cmdline_options(
      // this function is meant to require 3 fields: Phi, T, Conc
      const std::vector<string>& args,
      double& dt,
      int& Nt,
      int& Nv,
      double& hh_x,
      double& ww,
      double& shape_constant,
      double& mobility,
      double& kappa,
      double& c_alpha,
      double& c_beta,
      int& write_period,
      int_flags& flags,
      string& output_prefix,
      string& input_field_name,
      string& datasetPathPhi,
      string& datasetPathT,
      string& datasetPathConc,
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
      std::string datasetPath;
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
                     hh_x,
                     ww,
                     shape_constant,
                     mobility,
                     kappa,
                     c_alpha,
                     c_beta,
                     write_period,
                     output_prefix,
                     input_field_name,
                     datasetPath,
                     datasetPathPhi,
                     datasetPathT,
                     datasetPathConc,
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
      else if ( args[idx] == "-datasetPathPhi" )
      {
         if ( idx + 1 < args.size()) 
         {
            datasetPathPhi = string( args[idx + 1] );
            flags.datasetPathPhi = 1;  // true
            idx += 1;
         }
      }
      else if ( args[idx] == "-datasetPathT" )
      {
         if ( idx + 1 < args.size()) 
         {
            datasetPathT = string( args[idx + 1] );
            flags.datasetPathT = 1;  // true
            idx += 1;
         }
      }
      else if ( args[idx] == "-datasetPathConc" )
      {
         if ( idx + 1 < args.size()) 
         {
            datasetPathConc= string( args[idx + 1] );
            flags.datasetPathConc= 1;  // true
            idx += 1;
         }
      }
      else if ( args[idx] == "-datasetPath" )
      {
         if ( idx + 1 < args.size()) 
         {
            datasetPath = string( args[idx + 1] );
            flags.datasetPath = 1;  // true
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
      else if ( args[idx] == "-debug" ) 
      {
         flags.debug = 1;
      }
      else if ( args[idx] == "-mesh-size" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> hh_x;
            idx += 1;
         }
      }
      else if ( args[idx] == "-order-energy" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> ww;
            idx += 1;
         }
      }
      else if ( args[idx] == "-shape-constant" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> shape_constant;
            idx += 1;
         }
      }
      else if ( args[idx] == "-mobility" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> mobility;
            idx += 1;
         }
      }
      else if ( args[idx] == "-kappa" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> kappa;
            idx += 1;
         }
      }
      else if ( args[idx] == "-c-alpha" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> c_alpha;
            idx += 1;
         }
      }
      else if ( args[idx] == "-c-beta" )
      {
         if ( idx + 1 < args.size()) 
         {
            istringstream( args[idx + 1] ) >> c_beta;
            idx += 1;
         }
      }
   }

   if ((flags.datasetPathPhi == 0) || (flags.datasetPathT == 0)
            || (flags.datasetPathConc == 0))

   {
      if ( mynode == rootnode )
         cout << "Error: insufficient arguments;" 
            << " require dataset paths for three fields Phi, T, Conc"
            << " (-datasetPathPhi </phi>,"
           << " -datasetPathT </T>, -datasetPathConc </conc>)"
            << " and output prefix (-o <path/file> )" << endl;
      return EXIT_FAILURE;
   }
   else 
   {
      flags.datasetPath = 1; // in case it's needed by subsequent functions
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


int SPF_NS::read_cmdline_options_for_escape_time_ensembles(
      const std::vector<string>& args,
      double& dt,
      int& Nt,
      size_t& Nensemble,
      double& rate_scale_factor,
      int_flags& flags,
      string& output_prefix,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   string parameter_filename;

   // dummy variables used only for calling functions that require them
   int write_period,Nv;
   string input_field_name;

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
            //if ( mynode == rootnode )
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
      else if ( args[idx] == "-debug" ) 
      {
         flags.debug = 1;
      }
   }

   if ( !((flags.output_prefix == 1) & (flags.input_field == 1) ))
   {
      //if ( mynode == rootnode )
         cout << "Error: insufficient arguments;" 
            << " require both input field (-i <file>)"
            << " and output prefix (-o <path/file> )" << endl;

      return EXIT_FAILURE;
   }
   if ( flags.dt == 0 )
   {
      //if ( mynode == rootnode )
         cout << "warning: time increment (-dt <#>) not specified,"
           << " assigning dt = 1"
            << endl;
      dt = 1.0;
   }
   if ( flags.Nt == 0 )
   {
      //if ( mynode == rootnode )
         cout << "Error: number of time steps (-Nt <#> ) not specified"
            << endl;
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
} // read_cmdline_options_for_escape_time_ensembles()
#endif
