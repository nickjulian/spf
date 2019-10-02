/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: read_cmdline.hpp

#ifndef READ_CMDLINE_HPP
#define READ_CMDLINE_HPP

#define PRINT_USAGE cout << "Usage: " << argv[0] << " <OPTIONS>" << endl << "OPTIONS : " << endl << "   -i <input field hdf5 file>" << endl << "   -o <output file prefix>" << endl << "   -Nt <number of time steps>" << endl << "   [-dt <time increment>]" << endl;

#include <mpi.h>
#include <iostream>  // cout, cin, cerr, endl
#include <iomanip>   // setw, setprecision
//#include <fstream>   // ifstream, ofstream
#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <vector>
#include <string>

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
      string& output_prefix,
      string& input_field_name,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      );

} // SPF_NS
#endif
