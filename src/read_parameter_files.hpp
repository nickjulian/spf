/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: read_parameter_files.hpp

#ifndef READ_PARAMETER_FILES_HPP
#define READ_PARAMETER_FILES_HPP

#include <iostream>  // cout, cin, cerr, endl
#include <iomanip>   // setw, setprecision
#include <fstream>   // ifstream, ofstream
#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <string>
#include <vector>

using std::cout;
using std::cerr;

namespace SPF_NS
{
   //int read_parameter_file(
   //      );

   int read_state_file_list(
         const std::string& list_file_name,
         std::vector<std::string>& state_file_names
         );

}  // SPF_NS

#endif
