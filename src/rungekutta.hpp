/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: rungekutta.hpp
// Purpose:

#ifndef RUNGEKUTTA_HPP
#define RUNGEKUTTA_HPP

//#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE

namespace TEM_NS
{
   double marcusRK4( double (*integrand)(const double& tt, const double& yy,
                                    const double& poissonJump), 
                     const double& dt, 
                     const double& y0, 
                     const double& t0,
                     const double& jump);
}

#endif
