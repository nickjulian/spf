/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
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
