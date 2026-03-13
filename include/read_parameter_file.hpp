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
// File: read_parameter_file.hpp

#ifndef READ_PARAMETER_FILE_HPP
#define READ_PARAMETER_FILE_HPP

#define PRINT_USAGE std::cout << "Error, missing the parameter_file." << std::endl  \
                              << "Usage: mpirun -n #cores " << argv[0] << std::endl \
                              << "OPTIONS IN parameter_file: " << std::endl \
                              << "   -i <input field hdf5 file>" << std::endl \
                              << "   -o <output file prefix>" << std::endl \
                              << "   -datasetPath <concentration path to initial field within hdf5 file>" << std::endl \
                              << "   -TPath <temperature path to initial field within hdf5 file>" << std::endl \
                              << "   -fix_y <periodic or non-periodic boundary condition along y direction>" << std::endl \
                              << "   -c_extreme  <mixture contraction free energy adjustment for temperature, only use it when the final product has no phase separation>" << std::endl \
                              << "   -Nt <number of time steps>" << std::endl \
                              << "   -dt <time increment>" << std::endl \
                              << "   -wp <steps between file writes>" << std::endl \
                              << "   -c_fluc_theta <concentration fluctuation coefficient nu>" << std::endl \
                              << "   -mobility <mobility M>" << std::endl \
                              << "   -mesh_size <one voxel length (assuming cubic)>" << std::endl \
                              << "   -shape_constant <shape constant gamma>" << std::endl \
                              << "   -c_alpha <double-well free energy two well depth>" << std::endl \
                              << "   -c_beta <double-well free energy two well depth>" << std::endl \
                              << "   -kappa <gradient energy coefficient>" << std::endl \
                              << "   -h_d <thermal diffusivity>" << std::endl \
                              << "   -h_c <heat capacity>" << std::endl \
                              << "   -molar_volume <molar volume>" << std::endl \
                              << "   -Rho_W <W mass density>" << std::endl \
                              << "   -Rho_Cr <Cr mass density>" << std::endl \
                              << "   [-stat] <for debug>"<< std::endl \
                              << "   [-debug] <for debug>"<< std::endl;
//<< "   [-f <scalar integrand>]" << endl 
//<< "   [-r <rate scale factor>]" << endl 

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
         double& hh_x,
         double& shape_constant,
         double& mobility,
         double& kappa,
         double& c_alpha,
         double& c_beta,
         double& h_d,
         double& h_c,
         int& fix_y,
         int& c_extreme,
         double& c_fluc_theta,
         double& molar_volume,
         int& Rho_W,
         int& Rho_Cr,
         int& write_period,
         string& output_prefix,
         string& input_field_name,
         string& datasetPath,
         string& TPath,
         const int& mynode,
         const int& rootnode,
         MPI_Comm comm
         );

}  // SPF_NS

#endif
