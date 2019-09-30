/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: check_for_failure.cpp

#ifndef CHECK_FOR_FAILURE_CPP
#define CHECK_FOR_FAILURE_CPP
#include "check_for_failure.hpp"

bool SPF_NS::check_for_failure(
      int failflag,
      MPI_Comm comm
      )
{
   MPI_Allreduce( 
         MPI_IN_PLACE, // sendbuf
         &failflag, // recvbuf
         1, // count
         MPI_INT, // MPI datatype
         MPI_SUM, // reduction operation
         comm);
   if ( failflag != 0) return true;
   else return false;
}
#endif
