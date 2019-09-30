/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: jumps.cpp
// Purpose:

#ifndef JUMPS_CPP
#define JUMPS_CPP

#include "jumps.hpp"

#include <iostream>
using std::cerr;
using std::endl;


double TEM_NS::jumpDestination( const double& rr, const double& cc)
{
   return cc * exp(rr);
}

#endif
