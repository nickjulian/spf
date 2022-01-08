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
// File: stochastic_rates.hpp

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
         //const double& rate_scale_factor,
         const size_t& idx
         );

   double double_well_srscd(
            const double& c_t,
            const double& c_alpha,
            const double& c_beta,
            const double& lap, // laplacian
            const double& shape_constant,
            const double& gradient_coefficient
         );

   double double_well_mu_finel(
         const double& xx, // current voxel concentration
         const double& yy, // a neighbor's concentration
         const double& upward_shift,
         const double& ww,
         const double& boltz,
         const double& AA,
         const double& BB
         );

   double double_well_muhomo_finel(
         const double& xx,  
         const double& upward_shift,
         const double& ww,
         const double& boltz
         );

   int simple_identity_rate_gradient(
         double& local_rate,  // output
         const std::vector<double>& local_field,
         const size_t& neigh_idx,
         const size_t& idx
         );

   int simple_identity_rate_derivative( // (1/6)*rate_scale_factor
         double& local_rate_derivative,
         const std::vector<double>& local_field,
         //const double& rate_scale_factor,
         const size_t& idx
         );

   int simple_identity_rate_gradient_derivative( 
         double& local_rate_derivative,
         const std::vector<double>& local_field,
         const size_t& neigh_idx,
         const size_t& idx
         );

   double double_well_tilted(
         //double& local_rate,  // output
         //const std::vector<double>& local_field,
         const double& xx, // local_field[idx],
         //const double& rate_scale_factor,
         const double& upward_shift,
         const double& ww,
         const double& TT,
         const double& alpha
         //const size_t& idx
         );

   double double_well_tilted_gradient(
         const double& xx, // local_field[idx],
         const double& yy, // local_field[neigh_idxs[mm]],
         const double& upward_shift,
         const double& ww,
         const double& TT,
         const double& alpha
         );

   double double_well_tilted_derivative(
         //double& local_rate_derivative,  // output
         //const std::vector<double>& local_field,
         const double& xx, // local_field[idx],
         //const double& rate_scale_factor,
         const double& upward_shift,
         const double& ww,
         const double& TT,
         const double& alpha
         //const size_t& idx
         );

   double double_well_tilted_gradient_derivative(
         const double& xx, // local_field[idx],
         const double& yy, // local_field[idx],
         const double& upward_shift,
         const double& ww,
         const double& TT,
         const double& alpha
         );
} // SPF_NS
#endif
