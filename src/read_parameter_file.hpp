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
// File: read_parameter_file.hpp

#ifndef READ_PARAMETER_FILE_HPP
#define READ_PARAMETER_FILE_HPP

#include <mpi.h>
#include <iostream>  // cout, cin, cerr, endl
#include <iomanip>   // setw, setprecision
#include <fstream>   // ifstream, ofstream
#include <sstream>   // ostringstream
#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <string>
#include <vector>
#include <algorithm> // transform
#include <cctype>    // tolower

#include "flags.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::istringstream;
using std::ifstream;

namespace SPF_NS
{
   int read_parameter_file(
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
         );

   int read_parameter_file(
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
         );

   int read_parameter_file(
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
         );

   int read_parameter_file(
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
         string& input_field_file_name,
         string& datasetPath,
         string& datasetPathPhi,
         string& datasetPathT,
         string& datasetPathConc,
         const int& mynode,
         const int& rootnode,
         MPI_Comm comm
         );

   int read_parameter_file(
         const string& parameter_filename,
         int_flags& flags,
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
         string& output_prefix,
         string& input_field_file_name,
         string& datasetPath,
         string& datasetPathPhi,
         string& datasetPathT,
         string& datasetPathConc,
         const int& mynode,
         const int& rootnode,
         MPI_Comm comm
         );

   int read_state_file_list(
         const std::string& list_file_name,
         std::vector<std::string>& state_file_names
         );

}  // SPF_NS

#endif
