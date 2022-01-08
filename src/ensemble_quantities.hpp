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
// File: ensemble_quantities.hpp

#ifndef ENSEMBLE_QUANTITIES_HPP
#define ENSEMBLE_QUANTITIES_HPP

#include <math.h>    // log
//#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE

namespace SPF_NS
{
   int entropy_boltzmann( 
                     const int& Nv,
                     const int& population,
                     double& entropy);  // output
}

#endif
