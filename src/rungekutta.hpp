/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: rungekutta.hpp
// Purpose:

#ifndef RUNGEKUTTA_HPP
#define RUNGEKUTTA_HPP

//#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE

namespace SPF_NS
{
   double marcusRK4( double (*integrand)(const double& tt, 
                                          const double& yy,
                                          const double& poissonJump), 
                     const double& dt, 
                     const double& y0, 
                     const double& t0,
                     const double& jump);

   int marcusIntegral( double (*marcusIntegrand)(const double& tt,
                                          const double& yy,
                                          const double& dP), 
                     const double& rk_dt,    // increment of unit interval
                     const double& y0,    // path.back()
                     const double& t0,    // t_n, t_{n+1}=t_n + dt
                     const double& dP, // Y_{k} size of jump before marcus
                     double& jumpDestination);
}

#endif
