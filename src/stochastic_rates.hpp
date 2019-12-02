/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: stochastic_rates.hpp
// Purpose:

#ifndef STOCHASTIC_RATES_HPP
#define STOCHASTIC_RATES_HPP

#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <vector>
#include <iostream>  // cout, cin, cerr, endl
#include <algorithm> // sort
#include <math.h>    // sqrt

namespace SPF_NS
{
   int simple_identity_rate(// (1/6)*rate_scale_factor * local_field[idx]
         double& local_rate,
         const std::vector<double>& local_field,
         const double& rate_scale_factor,
         const size_t& idx
         );

   int simple_identity_rate_derivative( // (1/6)*rate_scale_factor
         double& local_rate_derivative,
         const std::vector<double>& local_field,
         const double& rate_scale_factor,
         const size_t& idx
         );

   int double_well_tilted(
         double& local_rate,  // output
         //const std::vector<double>& local_field,
         const double& xx, // local_field[idx],
         const double& rate_scale_factor,
         const double& ww,
         const double& TT,
         const double& alpha
         //const size_t& idx
         );

   int double_well_tilted_derivative(
         double& local_rate_derivative,  // output
         //const std::vector<double>& local_field,
         const double& xx, // local_field[idx],
         const double& rate_scale_factor,
         const double& ww,
         const double& TT,
         const double& alpha
         //const size_t& idx
         );
} // SPF_NS
#endif
