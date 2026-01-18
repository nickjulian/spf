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
// File: rand.hpp

#ifndef RAND_HPP
#define RAND_HPP

#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE
#include <random>    // default_random_engine, uniform_real_distribution
#include "Random123/philox.h"   // philox random number generator
#include "Random123/uniform.hpp"  // uniform distribution

#ifndef PI
#define PI 3.141592653589793238462643382795028814971693993751058209
#endif
#ifndef KB
#define KB 1.380649E-23
#endif

namespace SPF_NS
{
   class random
   {
      public:
         std::random_device rd;
         std::mt19937 generator;
         std::uniform_real_distribution<double> uniform_scale;

         using philox_type = r123::Philox4x64;
         philox_type rng;
         philox_type::key_type key;
         philox_type::ctr_type ctr;

         uint64_t counter = 0;

         // CONSTRUCTOR
         random(uint64_t seed = 0, uint64_t rank = 0)
         {
            generator = std::mt19937(rd());
            uniform_scale = std::uniform_real_distribution<double>( 0.0, 1.0);

            key[0] = seed;
            key[1] = rank;
            ctr = {{0, 0, 0, 0}};
         }

         // Reset Philox counter for a new timestep
         void reset_counter(uint64_t timestep)
         {
            counter = 0;
            ctr[1] = timestep; // vary RNG per timestep
         }

         // Generate uniform (0,1] using Philox
         uint64_t philox_generator()
         {
            ctr[0] = counter++; // bump counter
            auto result = rng(ctr, key);
            return result[0];
         }
   };
}
#endif
