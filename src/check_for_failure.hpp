/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: check_for_failure.hpp
// Purpose:

#ifndef CHECK_FOR_FAILURE_HPP
#define CHECK_FOR_FAILURE_HPP

#include <mpi.h>
#include "flags.hpp"

namespace SPF_NS 
{
   bool check_for_failure(
         int_flags& flags,
         MPI_Comm comm
         );
} // SPF_NS
#endif
