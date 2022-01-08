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
// File: read_cmdline.hpp

#ifndef READ_CMDLINE_HPP
#define READ_CMDLINE_HPP

#define PRINT_USAGE cout << "Usage: " << argv[0] << " <OPTIONS>" << endl << "OPTIONS : " << endl << "   -i <input field hdf5 file>" << endl << "   -datasetPath <path to initial field within hdf5 file>" << endl << "   -o <output file prefix>" << endl << "   -Nt <number of time steps>" << endl << "   -Nv <max number of states per voxel>" << endl << "   [-dt <time increment>]" << endl << "   [-wp <steps between file writes>]" << endl << "   [-stat]"<< endl;
//<< "   [-f <scalar integrand>]" << endl 
//<< "   [-r <rate scale factor>]" << endl 

#include <mpi.h>
#include <iostream>  // cout, cin, cerr, endl
#include <iomanip>   // setw, setprecision
//#include <fstream>   // ifstream, ofstream
#include <sstream>   // ostringstream
#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <vector>
#include <string>

#include "flags.hpp" // int_flags
#include "read_parameter_file.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::setw;
using namespace std;

namespace SPF_NS
{

int read_cmdline_options(
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
      );

int read_cmdline_options(
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
      );


int read_cmdline_options(
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
      );

int read_cmdline_options(
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
      );

int read_cmdline_options(
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
      );

int read_cmdline_options(
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
      );

int read_cmdline_options(
      const std::vector<string>& args,
      double& dt,
      int& Nt,
      int& write_period,
      int_flags& flags,
      string& output_prefix,
      string& input_field_name,
      string& datasetPath
      );

int read_cmdline_options(
      const std::vector<string>& args,
      double& dt,
      int& Nt,
      int& write_period,
      int_flags& flags,
      string& output_prefix,
      string& input_field_name
      );

int read_cmdline_options_for_escape_time_ensembles(
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
      );
} // SPF_NS
#endif
