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
// File: read_parameter_file.cpp

#ifndef READ_PARAMETER_FILE_CPP
#define READ_PARAMETER_FILE_CPP

#include "read_parameter_file.hpp"

int SPF_NS::read_parameter_file(
    const string &parameter_filename,
    int_flags &flags,
    double &dt,
    int &Nt,
    double &hh_x,
    double &shape_constant,
    double &mobility,
    double &kappa,
    double &c_alpha,
    double &c_beta,
    double &h_d,
    double &h_c,
    int &fix_y,
    int &c_extreme,
    double &c_fluc_theta,
    double &molar_volume,
    int &Rho_W,
    int &Rho_Cr,
    int &write_period,
    string &output_prefix,
    string &input_field_name,
    string &datasetPath,
    string &TPath,
    const int &mynode,
    const int &rootnode,
    MPI_Comm comm)
{
   ifstream parameter_file( parameter_filename.c_str() );
   //flags.parameter_file = 1;
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
         //else if (! line_chunk.compare("-Nv") )
         //{
         //   file_line_stream >> Nv;
         //   flags.Nv = 1;
         //}
         else if (! line_chunk.compare("-mesh_size") )
         {
            file_line_stream >> hh_x;
         }
         //else if (! line_chunk.compare("-order-energy") )
         //{
         //   file_line_stream >> ww;
         //}
         else if (! line_chunk.compare("-shape_constant") )
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
         else if (! line_chunk.compare("-c_alpha") )
         {
            file_line_stream >> c_alpha;
         }
         else if (! line_chunk.compare("-c_beta") )
         {
            file_line_stream >> c_beta;
         }
         else if (! line_chunk.compare("-h_d") )
         {
            file_line_stream >> h_d;
         }
         else if (! line_chunk.compare("-h_c") )
         {
            file_line_stream >> h_c;
         }
         else if (! line_chunk.compare("-wp") )
         {
            file_line_stream >> write_period;
            flags.wp = 1;
         }
         else if (! line_chunk.compare("-fix_y") )           
         {
            file_line_stream >> fix_y;
         }
         else if (!line_chunk.compare("-c_extreme"))
         {
            file_line_stream >> c_extreme;
         }
         else if (! line_chunk.compare("-c_fluc_theta") )           
         {
            file_line_stream >> c_fluc_theta;
         }
         else if (! line_chunk.compare("-molar_volume") )           
         {
            file_line_stream >> molar_volume;
         }
         else if (! line_chunk.compare("-Rho_W") )           
         {
            file_line_stream >> Rho_W;
         }
         else if (! line_chunk.compare("-Rho_Cr") )           
         {
            file_line_stream >> Rho_Cr;
         }
         else if (!line_chunk.compare("-TPath"))               
         {
            file_line_stream >> TPath;
            flags.datasetPathT = 1;
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
      return EXIT_SUCCESS;
   }
   return EXIT_FAILURE;
}
#endif
