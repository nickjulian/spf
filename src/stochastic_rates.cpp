/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: stochastic_rates.cpp
// Purpose:

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
                  * sqrt(local_field[idx]);
   return EXIT_SUCCESS;
}

int SPF_NS::simple_identity_rate_derivative( 
      double& local_rate_derivative,
      const std::vector<double>& local_field,
      //const double& rate_scale_factor,
      const size_t& idx
      )
{
   if (local_field[idx] <= 0.0 ) 
         local_rate_derivative = 100000000000000000000.0; 
   //if (xx >= 1.0 ) return ONESIXTH * rate_scale_factor * 0.5;
   local_rate_derivative = ONESIXTH //* rate_scale_factor 
                              * 0.5/sqrt(local_field[idx]);
   return EXIT_SUCCESS;
}

double SPF_NS::double_well_tilted(
      //double* local_rate,  // output
      //const std::vector<double>& local_field,
      const double& xx, // local_field[idx],
      //const double& rate_scale_factor,
      const double& ww,
      const double& TT,
      const double& alpha
      //const size_t& idx
      )
{
   if (xx <= 0.0 ) return 0.0; // this is probably bad
   if (xx >= 1.0 ) return 1E20;  // huge number
   // \omega(2x-1)+k_{B}T[\alpha ln(x/(1-x))+(alpha-1)/(1-x)]
   return ww *(2*xx -1) + (0.00008617)*TT*(
               alpha * log(xx/(1.0 - xx)) + (alpha -1)/(1 - xx));
}

double SPF_NS::double_well_tilted_derivative(
      //double& local_rate_derivative, // output
      //const std::vector<double>& local_field,
      const double& xx, // local_field[idx],
      //const double& rate_scale_factor,
      const double& ww,
      const double& TT,
      const double& alpha
      //const size_t& idx
      )
{
   if (xx <= 0.0 ) return 0.0; // this is probably bad
   if (xx >= 1.0 ) return -1E20;  // huge negative number
   //double xx; xx = local_field[idx];
   // 2*\omega+k_{B}T[(\alpha/x) + (\alpha/(1-x)) - (\alpha -1)/((1-x)^{2})]
   //if ( (xx == 0.0 ) || (xx == 1.0))
   //   local_rate_derivative = 

   //if ( (xx <= 0.0 ) || (xx >= 1.0))
   //{
   //   return xx;
   //}
   double sigma;
   sigma = sqrt(ww *(2*xx -1) + (0.00008617)*TT*(
               alpha * log(xx/(1.0 - xx)) + (alpha -1)/(1 - xx)));

   return ((0.5/sigma)
      *(2.0*ww + (0.00008617)*TT*( 
                                    (alpha/xx) + (alpha/(1.0-xx)) 
                                    + (alpha -1.0)/((1.0-xx)*(1.0-xx))
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
