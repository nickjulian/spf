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
// File: stochastic_rates.cpp

#ifndef STOCHASTIC_RATES_CPP
#define STOCHASTIC_RATES_CPP

#include "stochastic_rates.hpp"

#ifndef ONESIXTH 
#define ONESIXTH 0.16666666666666666666666666666666666666666666666666666667
#endif

int SPF_NS::simple_identity_rate(  // rate_scale_factor * local_field[idx]
      double& local_rate,  // output
      const std::vector<double>& local_field,
      //const double& rate_scale_factor,
      const size_t& idx
      )
{
   if (local_field[idx] <= 0.0 ) local_rate = 0;
   local_rate = ONESIXTH //* rate_scale_factor 
                  * (local_field[idx]);
   return EXIT_SUCCESS;
}

int SPF_NS::simple_identity_rate_derivative( 
      double& local_rate_derivative,
      const std::vector<double>& local_field,
      //const double& rate_scale_factor,
      const size_t& idx
      )
{  // has the sqrt() incorporated because it's only used to drive Gaussian 
   if (local_field[idx] <= 0.0 ) 
         local_rate_derivative = 100000000000000000000.0; 
   //if (xx >= 1.0 ) return ONESIXTH * rate_scale_factor * 0.5;
   local_rate_derivative = ONESIXTH //* rate_scale_factor 
                              * 0.5/sqrt(local_field[idx]);
   return EXIT_SUCCESS;
}

int SPF_NS::simple_identity_rate_gradient(  // rate_scale_factor * local_field[idx]
      double& local_rate,  // output
      const std::vector<double>& local_field,
      const size_t& neigh_idx,
      const size_t& idx
      )
{
   //if ((local_field[idx] < 0.0 ) || (local_field[neigh_idx] < 0.0))
   //{
   //   local_rate = 0;
   //}
   //else
   //{
   local_rate = ONESIXTH // TODO: divide by voxel separation
                  * (
                        local_field[idx]
                        - local_field[neigh_idx]
                        );
   //}
   return EXIT_SUCCESS;
}

int SPF_NS::simple_identity_rate_gradient_derivative( 
      double& local_rate_derivative,
      const std::vector<double>& local_field,
      const size_t& neigh_idx,
      const size_t& idx
      )
{
   if( local_field[idx] > local_field[neigh_idx])
   {
      local_rate_derivative = ONESIXTH //* rate_scale_factor 
                                 * 0.5/sqrt(
                                       local_field[idx]
                                          - local_field[neigh_idx]
                                       );
   }
   else
   if( local_field[idx] < local_field[neigh_idx])
   {  // TODO: ensure this is correct
      local_rate_derivative = ONESIXTH //* rate_scale_factor 
                                 * -0.5/sqrt(
                                       local_field[neigh_idx]
                                       - local_field[idx]
                                       );
   }
   else  //  gradient == 0
   {
      local_rate_derivative = 100000000000000000000.0; 
   }
      
   return EXIT_SUCCESS;
}

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
            *(1 - 2*c_t) - gradient_coefficient*lap;
}

double SPF_NS::double_well_muhomo_finel(
      const double& xx,  
      const double& upward_shift,
      const double& ww,
      const double& boltz
      )
{
   // homogeneous part of alloy chemical potential
   if (xx <= 0.0 ) return 0.0; // this is probably bad
   if (xx >= 1.0 ) return 1E10;  // huge number
   return ww * (2*xx - 1.0)
            + boltz*log(xx/(1.0 - xx))
            + upward_shift;
}

double SPF_NS::double_well_mu_finel(
      const double& xx, // current voxel concentration
      const double& yy, // a neighbor's concentration
      const double& upward_shift,
      const double& ww,
      const double& boltz,
      const double& AA,
      const double& BB
      )
{
   // alloy chemical potential including concentration gradient
   //  and lambda (coefficient of the gradient).
   return double_well_muhomo_finel(xx, upward_shift, ww, boltz)
            + (AA + (BB * xx * (1 - xx))) // lambda
               * abs(xx - yy);   // TODO: questionable gradient
}

double SPF_NS::double_well_tilted(
      //double* local_rate,  // output
      //const std::vector<double>& local_field,
      const double& xx, // local_field[idx],
      //const double& rate_scale_factor,
      const double& upward_shift,
      const double& ww,
      const double& TT,
      const double& alpha
      //const size_t& idx
      )
{
   if (xx <= 0.0 ) return 0.0; // this is probably bad
   if (xx >= 1.0 ) return 1E10;  // huge number
   // \omega(2x-1)+k_{B}T[\alpha ln(x/(1-x))+(alpha-1)/(1-x)]
   return ww *(2*xx -1) + 0.00008617*TT*(
                  alpha * log(xx/(1.0 - xx)) + (alpha -1)/(1 - xx)
               ) + upward_shift;
}

double SPF_NS::double_well_tilted_gradient(
      //double* local_rate,  // output
      //const std::vector<double>& local_field,
      const double& xx, // local_field[idx],
      const double& yy, // local_field[neigh_idxs[mm]],
      //const double& rate_scale_factor,
      const double& upward_shift,
      const double& ww,
      const double& TT,
      const double& alpha
      //const size_t& idx
      )
{
   double sgn; sgn = 1.0;
   if (xx > yy) sgn = 1.0;
   if (xx < yy) sgn = -1.0;
   if ((xx - yy) == 0.0 ) return 0.0; // this is probably bad
   if ((xx - yy) == 1.0 ) return 1E10;  // huge number
   if ((xx - yy) == -1.0 ) return -1E10;  // huge number
   // \omega(2x-1)+k_{B}T[\alpha ln(x/(1-x))+(alpha-1)/(1-x)]
   return 
      sgn*(
         ww *(2*sgn*(xx-yy) -1) + 0.00008617*TT*(
                  alpha * log(sgn*(xx-yy)/(1.0 - sgn*(xx-yy))) 
                  + (alpha -1)/(1 - sgn*(xx-yy))
               ) + upward_shift
         );
}

double SPF_NS::double_well_tilted_derivative(
      //double& local_rate_derivative, // output
      //const std::vector<double>& local_field,
      const double& xx, // local_field[idx],
      //const double& rate_scale_factor,
      const double& upward_shift,
      const double& ww,
      const double& TT,
      const double& alpha
      //const size_t& idx
      )
{
   if (xx <= 0.0 ) return 0.0; // this is probably bad
   if (xx >= 1.0 ) return 1E10;  // huge number
   //double xx; xx = local_field[idx];
   // 2*\omega+k_{B}T[(\alpha/x) + (\alpha/(1-x)) - (\alpha -1)/((1-x)^{2})]
   //if ( (xx == 0.0 ) || (xx == 1.0))
   //   local_rate_derivative = 

   //if ( (xx <= 0.0 ) || (xx >= 1.0))
   //{
   //   return xx;
   //}
   double sigma;
   double tmpdbl;
   tmpdbl = ww *(2*xx -1) + (0.00008617)*TT*(
                  alpha * log(xx/(1.0 - xx)) + (alpha -1)/(1 - xx)
               )
               + upward_shift;
   if ( tmpdbl > 0 ) sigma = sqrt( tmpdbl );
   else if (tmpdbl <= 0.0)
   { // chemical potential is < 0; the rate has been set to 0 elsewhere
      return 0.0;
   }

   return ((0.5/sigma)
      *(2.0*ww + (0.00008617)*TT*( 
                                    (alpha/xx) + (alpha/(1.0-xx))
                                    + (alpha -1.0)/((1.0-xx)*(1.0-xx))
                                 )));

}

double SPF_NS::double_well_tilted_gradient_derivative(
      //double& local_rate_derivative, // output
      //const std::vector<double>& local_field,
      const double& xx, // local_field[idx],
      const double& yy, // local_field[idx],
      //const double& rate_scale_factor,
      const double& upward_shift,
      const double& ww,
      const double& TT,
      const double& alpha
      //const size_t& idx
      )
{
   double sgn; sgn = 1.0;
   if (xx > yy) sgn = 1.0;
   if (xx < yy) sgn = -1.0;
   if ((xx - yy) == 0.0 ) return 0.0; // this is probably bad
   if ((xx - yy) == -1.0 ) return -1E10;  // huge number
   if ((xx - yy) == 1.0 ) return 1E10;  // huge number

   double sigma;
   double tmpdbl;
   tmpdbl = ww *(2*sgn*(xx-yy) -1) + (0.00008617)*TT*(
                  alpha * log(sgn*(xx - yy)/(1.0 - sgn*(xx-yy))) 
                  + (alpha -1)/(1 - sgn*(xx-yy))
               )
               + upward_shift;
   if ( tmpdbl > 0 ) sigma = sqrt( tmpdbl );
   else if (tmpdbl < 0.0)
   {
      sigma = sqrt( -1.0* tmpdbl);
   }
   else if ((tmpdbl == 0.0) || isnan(tmpdbl))
   {
      return 0.0;
   }

   return sgn*((0.5/sigma)
      *(2.0*ww + (0.00008617)*TT*( 
                                    (alpha/(sgn*(xx - yy))) 
                                    + (alpha/(1.0-sgn*(xx - yy))) 
                                    + (alpha -1.0)/((1.0-sgn*(xx - yy))
                                       *(1.0-sgn*(xx - yy)))
                                 )));
}

//int SPF_NS::double_well_tilted_ie(  // includes interface energy?
//      double& local_rate,  // output
//      const std::vector<double>& local_field,
//      const std::vector<size_t>& neigh_idxs,  // 6 elements
//      const double& rate_scale_factor,
//      const double& ww,
//      const double& TT,
//      const double& alpha,
//      const size_t& idx
//      )
//{
//   // \omega(2x-1)+k_{B}T[\alpha ln(x/(1-x))+(alpha-1)/(1-x)]  
//   //    + \beta \nabla x \nabla^{2} x
//   local_rate = ww *(2* local_field[idx] -1) + (8.617E-5)*TT*(
//         alpha * log(local_field[idx]/(1.0 - local_field[idx])) 
//         + (alpha -1)/(1 - local_field[idx]));
//   // TODO: add the gradient and laplacian ...
//   return EXIT_SUCCESS;
//}
#endif
