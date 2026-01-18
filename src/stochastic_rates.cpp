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
// File: stochastic_rates.cpp

#ifndef STOCHASTIC_RATES_CPP
#define STOCHASTIC_RATES_CPP

#include "stochastic_rates.hpp"

#ifndef ONESIXTH 
#define ONESIXTH 0.16666666666666666666666666666666666666666666666666666667
#endif

double SPF_NS::double_well_srscd(
         const double& c_t,
         const double& c_alpha,
         const double& c_beta,
         const double& lap, // laplacian
         const double& shape_constant,
         const double& gradient_coefficient
      )
{
   return 2*shape_constant*(c_t - c_alpha)*(c_beta - c_t)
            *(1 - 2*c_t) - gradient_coefficient*gradient_coefficient*lap;
}
#endif
