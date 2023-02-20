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
// File: flags.hpp

#ifndef FLAGS_HPP
#define FLAGS_HPP

#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE

namespace SPF_NS
{
struct int_flags
{
   // 0 == false;
   unsigned int dt;
   unsigned int Nt;
   unsigned int Nv;
   unsigned int wp;
   unsigned int datasetPath;
   unsigned int datasetPathPhi;
   unsigned int datasetPathT;
   unsigned int datasetPathConc;
   unsigned int output_prefix;
   unsigned int input_field;
   unsigned int calcstat;
   unsigned int parameter_file;
   unsigned int debug;
   unsigned int fail;

   int_flags() // constructor
   {
      dt = 0;  // 0 == false
      Nt = 0;
      Nv = 0;
      wp = 0;
      output_prefix = 0;
      input_field = 0;
      calcstat = 0;
      parameter_file = 0;
      debug = 0;
      fail = 0;
      datasetPath = 0;
      datasetPathPhi = 0;
      datasetPathT = 0;
      datasetPathConc = 0;
   }
};
} // SPF_NS

#endif
