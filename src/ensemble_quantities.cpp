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
// File: ensemble_quantities.cpp
// Purpose: implementation of various Runge-Kutta methods

#ifndef ENSEMBLE_QUANTITIES_CPP
#define ENSEMBLE_QUANTITIES_CPP

#include "ensemble_quantities.hpp"

#include <iostream> // cout, endl;

int SPF_NS::entropy_boltzmann( 
                     const int& Nv,
                     const int& population,
                     double& entropy)  // output
{
   // S = log(\Omega), \Omega = number of equivalent microstates per this 
   //    macrostate. Omitting Boltzmann's constant for now.
   
   // \Omega = Nv choose population
   int omega; omega = 1;
   size_t ii;
   if ( 2* population > Nv )
   {
      for ( ii=Nv; ii > population; --ii)
      {
         omega *= ii;
      }
      for ( ii= Nv - population ; ii > 1; --ii)
      {
         omega /= ii;
      }
   }
   else
   {
      for ( ii=Nv; ii > Nv - population; --ii)
      {
         omega *= ii;
      }
      for ( ii= population ; ii > 1; --ii)
      {
         omega /= ii;
      }
   }
   if ( omega <= 0)
   {
      std::cout << "Error, entropy_boltzmann, omega <= 0" << std::endl;
      return EXIT_FAILURE;
   }
   entropy = log( double(omega) );

   return EXIT_SUCCESS;
}


#endif
