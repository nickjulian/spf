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
// File: macheps.hpp

#ifndef MACHEPS_HPP
#define MACHEPS_HPP

namespace SPF_NS
{
class epsilon
{
   public:
   double dbl;
   double dblsqrt;
   //float flt;

   // constructor
   epsilon()
   {
      double eps; eps = 1.0;
      double trial; trial = 1.0 + eps;
      while( 1.0 != trial )
      {
         eps *= 0.5;
         trial = 1.0 + eps;
      }
      dbl = eps; 
      dblsqrt = sqrt(eps);
   }
};
} // SPF_NS
#endif
