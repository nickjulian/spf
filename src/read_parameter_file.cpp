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
// File: read_parameter_file.cpp

#ifndef READ_PARAMETER_FILE_CPP
#define READ_PARAMETER_FILE_CPP

#include "read_parameter_file.hpp"

int SPF_NS::read_state_file_list(
               const std::string& list_file_name,
               std::vector<std::string>& state_file_names
               )
{
   std::ifstream list_file( list_file_name.c_str() );
   if ( list_file.is_open() )
   {
      std::string state_file_name;
      std::string list_file_line;
      size_t first_pos; 
      // iterate over the lines in the file
      while (std::getline( list_file, list_file_line ) && list_file.good())
      {
         // skipping lines of whitespace
         first_pos = list_file_line.find_first_not_of(" \t");
         while ( first_pos == std::string::npos )
         {
            std::getline( list_file, list_file_line );
            first_pos = list_file_line.find_first_not_of(" \t");
         }
         std::istringstream list_file_line_stream( list_file_line );
         list_file_line_stream >> state_file_name;

         state_file_names.push_back( state_file_name );
      }
      //if ( ! list_file.good() )
      //{
      //   std::cerr << "Error, state of opened file "
      //      << list_file_name
      //      << ".good() is false" << std::endl;
      //   list_file.close();
      //   std::cout << "state_file_names : ";
      //   for ( std::vector<std::string>::const_iterator itr = state_file_names.begin(); itr != state_file_names.end(); ++itr)
      //      std::cout << std::endl << "   " << *itr;
      //   std::cout << std::endl;
      //   return EXIT_FAILURE;
      //}
   }
   else
   {
      std::cerr << "Error opening file " << list_file_name << std::endl;
      list_file.close();
      return EXIT_FAILURE;
   }
   list_file.close();
   return EXIT_SUCCESS;
}

int SPF_NS::read_parameter_file(
      const string& parameter_filename,
      int_flags& flags,
      double& dt,
      int& Nt,
      int& Nv,
      int& write_period,
      string& output_prefix,
      string& input_field_name,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
      double hh_x;
      double ww;
      double shape_constant;
      double mobility;
      double kappa;
      double c_alpha;
      double c_beta;
      std::string datasetPath = "/phi";
      return read_parameter_file(
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
                  );
}

int SPF_NS::read_parameter_file(
      const string& parameter_filename,
      int_flags& flags,
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
      string& output_prefix,
      string& input_field_name,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
      std::string datasetPath = "/phi";
      return read_parameter_file(
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
                  );
}

int SPF_NS::read_parameter_file(
      const string& parameter_filename,
      int_flags& flags,
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
      string& output_prefix,
      string& input_field_name,
      string& datasetPath,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   ifstream parameter_file( parameter_filename.c_str() );
   if ( parameter_file.is_open() )
   {
      //cout << "reading " << parameter_filename.c_str() << endl;
      string file_line;
      while( getline( parameter_file, file_line) && parameter_file.good() )
      {
         // skip empty whitespace lines
         size_t first = file_line.find_first_not_of(" \t");
         while( first == std::string::npos)
         {
            getline( parameter_file, file_line);
            first = file_line.find_first_not_of(" \t");
         }

         istringstream file_line_stream( file_line );
         string line_chunk;
         file_line_stream >> line_chunk;
         // make ASCII lower case
         //std::transform( line_chunk.begin(), line_chunk.end(),
         //            line_chunk.begin(), (int(*)(int))tolower );
         
         // compare line_chunks to cmdline flags
         //std::cout << "reading lines" << std::endl;//debug
         if ( ! line_chunk.compare("-o") )
         {
            file_line_stream >> output_prefix;
            flags.output_prefix = 1;
            //std::cout << "output_prefix : " << output_prefix << std::endl;//debug
         }
         else if (! line_chunk.compare("-i") )
         {
            file_line_stream >> input_field_name;
            flags.input_field = 1;
         }
         else if (! line_chunk.compare("-datasetPath"))
         {
            file_line_stream >> datasetPath;
            flags.datasetPath = 1;
         }
         else if (! line_chunk.compare("-dt") )
         {
            file_line_stream >> dt;
            flags.dt = 1;
         }
         else if (! line_chunk.compare("-Nt") )
         {
            file_line_stream >> Nt;
            flags.Nt = 1;
         }
         else if (! line_chunk.compare("-Nv") )
         {
            file_line_stream >> Nv;
            flags.Nt = 1;
         }
         else if (! line_chunk.compare("-mesh-size") )
         {
            file_line_stream >> hh_x;
         }
         else if (! line_chunk.compare("-order-energy") )
         {
            file_line_stream >> ww;
         }
         else if (! line_chunk.compare("-shape-constant") )
         {
            file_line_stream >> shape_constant;
         }
         else if (! line_chunk.compare("-mobility") )
         {
            file_line_stream >> mobility;
         }
         else if (! line_chunk.compare("-kappa") )
         {
            file_line_stream >> kappa;
         }
         else if (! line_chunk.compare("-c-alpha") )
         {
            file_line_stream >> c_alpha;
         }
         else if (! line_chunk.compare("-c-beta") )
         {
            file_line_stream >> c_beta;
         }
         else if (! line_chunk.compare("-wp") )
         {
            file_line_stream >> write_period;
            flags.wp = 1;
         }
         else if (! line_chunk.compare("-stat") )
         {
            flags.calcstat = 1;
         }
         else if (! line_chunk.compare("-debug") )
         {
            flags.debug = 1;
         }
         else
         {
            if ( mynode == rootnode )
            {
               cerr << "Error, unexpectied argument in parameter file: "
                  << file_line << endl;
            }
            parameter_file.close();
            return EXIT_FAILURE;
         }
      }
   }
   return EXIT_SUCCESS;
}

#endif
