/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
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
   }
};
} // SPF_NS

#endif
