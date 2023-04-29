
// File: free_energy_gl.hpp
//
#ifndef FREE_ENERGY_GL_HPP
#define FREE_ENERGY_GL_HPP
#include <cstdlib>
#include <cassert>
#include <cmath> // erf() (maybe try <math.h> or <ctgmath>)
#include <vector>
#include <list>
//#include <iostream> // std::cout // debug

#include "rungekutta.hpp"

#ifndef ONESIXTH 
#define ONESIXTH 0.16666666666666666666666666666666666666666666666666666667
#endif

namespace SPF_NS
{
class free_energy_gl
{
   // Equations
   // 
   // F(conc,phi,T)
   // = c1(1-phi)^2 (conc-(1-c2(erf((T-c3)/c4)+1)))^2 
   //   + c1 phi^2 (c2(1+erf((T-c3)/c4)) - conc)^2 + c5 \phi^4 
   //   + c6*(\nabla phi)^2   // note that this is not \nabla^{2} phi
   // c1 = 100    = Fconstants[0]
   // c2 = 0.25   = Fconstants[1]
   // c3 = 7      = Fconstants[2]
   // c4 = 3      = Fconstants[3]
   // c5 = 0      = Fconstants[4]
   // c6 = 0      = Fconstants[5]
   // \nabla phi = [dphi/dx1, dphi/dx2, ...], is a vector
   // 
   // dT = D_{T(t)} \nabla^2 T(t) dt + (L/c_P(T(t))) \diamond dphi
   // dphi = mu(t) dt + \sqrt{mu(t) dt} \circle dB(t)
   // \mu = -M_{phi} dF/dphi
   // \Delta conc_{i\rightarrow j}(t) = N_{lambda_i}(\Delta t, 1)
   // with rate lambda_i = \M_c \nabla^2 dF/dconc
   // 
   // dF if only phi may vary ( to find dF/dphi)
   // = -2 c1(1-phi) (conc-(1-c2(erf((T-c3)/c4)+1)))^2 dphi
   //   + 2 c1 phi (c2(1+erf((T-c3)/c4)) - conc)^2 dphi
   //   + 4 c5 \phi^3 dphi + c6*d( \nabla^{2] phi)/dphi
   //
   // d/dphi (\nabla^{2} phi)
   //    \approx
   //       d/dphi (1/hh_x)((phi_{x-1} -2*phi_{x} +phi_{x+1}))
   //       d/dphi (1/hh_y)((phi_{y-1} -2*phi_{y} +phi_{y+1}))
   //       d/dphi (1/hh_z)((phi_{z-1} -2*phi_{z} +phi_{z+1}))
   //    \approx -2*(1/hh_x + 1/hh_y + 1/hh_z)     ?????????

   // dF/dconc
   // = 2 c1 (1-phi)^2 (conc-(1-c2(erf((T-c3)/c4)+1)))
   //   - 2 c1 phi^2 (c2(1+erf((T-c3)/c4)) - conc)

   private:
      double Fconstants[6]; // used in dF
      double prefactor; // cPrefactor, used in C_{P}(T)
      double alpha; // used in C_{P}(T)
      double T0; // used in Keq, C_{P}(T)
      double cbase; // used in C_{P}(T)
      double latentHeat; // used in dT
      double mobility_phi[3]; // used in dphi
      double mobility_conc[3]; // used in dconc
      double orderEnergy; // used in delta_phi
      std::vector<double> unitInterval;
      double rk_dt;
   public:
      double thermalDiffusivity; // used in dT, and stability check
      double kb;
      // CONSTRUCTOR
      free_energy_gl(
            const double& c1, // used in F(conc,phi,T)
            const double& c2, // used in F(conc,phi,T)
            const double& c3, // used in F(conc,phi,T)
            const double& c4, // used in F(conc,phi,T)
            const double& c5, // used in F(conc,phi,T)
            const double& c6, // used in F(conc,phi,T)
            const double& cPrefactor,  // used in C_{P}(T), -0.05
            const double& aa, // alpha, used in C_{P}(T)
            const double& TT, // T0, used in C_{P}(T), 10
            const double& cb, // cbase, used in C_{P}(T), 2
            const double& L, // latentHeat
            const double& ww, // order energy for phase
            const double& M_phi_x, // mobility of phase
            const double& M_phi_y, // 3 components required due to 
            const double& M_phi_z, // dependency on grid spacing
            const double& M_conc_x, // mobility of conc. (negl. seperation)
            const double& M_conc_y, // mobility of conc. (negl. seperation)
            const double& M_conc_z, // mobility of conc. (negl. seperation)
            const double& D_T,
            const double& rk_dtIn
            )
      {
         Fconstants[0] = c1;
         Fconstants[1] = c2;
         Fconstants[2] = c3;
         Fconstants[3] = c4;
         Fconstants[4] = c5;
         Fconstants[5] = c6;
         prefactor = cPrefactor;
         alpha = aa;
         T0 = TT;
         cbase = cb;
         latentHeat = L;
         orderEnergy = ww;
         mobility_phi[0] = M_phi_x;
         mobility_phi[1] = M_phi_y;
         mobility_phi[2] = M_phi_z;
         mobility_conc[0] = M_conc_x;
         mobility_conc[1] = M_conc_y;
         mobility_conc[2] = M_conc_z;
         thermalDiffusivity = D_T;
         rk_dt = rk_dtIn;
         kb = 1.0;
         //kb = 0.31320406192891;
         for ( size_t t=0; t < 1.0/rk_dt; ++t)
         {
            unitInterval.push_back(t * rk_dt);
         }
      }
      //// DESTRUCTOR
      //~free_energy_gl()
      //{
      //}
      // PUBLIC MEMBER FUNCTIONS
      double F( 
            const double& conc, // in units of # walkers
            const double& Nv,  // to convert to concentration
            const double& phi,
            const double& TT,
            const double& lap1dPhiX,
            const double& lap1dPhiY,
            const double& lap1dPhiZ
            );

      double dF_dconc(
            const double& conc, // in units of # of walkers
            const double& Nv, // for conversion to/from walkers/conc.
            const double& phi,
            const double& T
            );
      //{
      //   // dF/dconc
      //   // = 2 c1 (1-phi)^2 (conc-(1-c2(erf((T-c3)/c4)+1)))
      //   //   - 2 c1 phi^2 (c2(1+erf((T-c3)/c4)) - conc)
      //   return 
      //      2 * Fconstants[0] * (1-phi)*(1-phi)
      //       * (conc -(1-Fconstants[1]
      //                * (erf((T-(Fconstants[2]))/(Fconstants[3]))+1))
      //         )
      //      - 2 * Fconstants[0] * phi * phi
      //       * (Fconstants[1]
      //                * (1+erf((T-(Fconstants[2]))/(Fconstants[3])))
      //                - conc); // + 0 + 0
      //}

      double entropy( 
            const double& conc, // in units of # of walkers
            const double& phi
            );

      double jump_rate_conc(
                           const double& lap1d,
                           const double& hh,
                           const size_t& nn)
      {
         // -M_conc \nabla^2 dF/dconc
         //return (-1*mobility_conc[nn]/hh)*lap1d;// hh is in mobility_conc
         return (-1*mobility_conc[nn])*lap1d;
         // shouldn't order energy be here too?
      }

      double delta_T_decoupled(
            //  \frac{\partial T}{\partial t} = D_{T} \nabla^{2} T
            const double& TT,
            const double& lapT, // laplacian of T
            const double& deltaPhi,
            const double& dt
            );

      double delta_T(
            //  \frac{\partial T}{\partial t} = D_{T} \nabla^{2} T
            const double& TT,
            const double& lapT, // laplacian of T
            const double& deltaPhi,
            const double& dt
            );

      double delta_T_marcus(
            //double (*latent_heat_c_p_inverse)(
            //   const double& tt,
            //   const double& TT,
            //   const double& jump
            //   ),
            const double& TT,
            const double& lapT,
            const double& dt,
            const double& phi_initial,
            const std::list<double>& deltaPhiSequence //sequence of jumps in Phi
            );

      int delta_phi_marcus(
            std::list<double>& deltaPhiSequence, // output
            const int& conc_local_flux_sum,
            const double& conc_initial,
            const double& Nv, // for conversion to/from walkers/conc.
            const double& phi,
            const double& lap3dPhi,
            const double& TT,
            const std::vector<double>& hh
            );

      double delta_phi_noiseless(
            // -M_{phi} \frac{\delta F}{\delta \phi}
            const double& concs, // in units of # of walkers
            const double& Nv, // for conversion to/from walkers/conc.
            const double& phi,
            const double& lap3dPhi,
            const double& TT,
            const std::vector<double>& hh
            );

      double marcusRK4(
            const double& dt,
            const double& y0,
            const double& t0,
            const double& jump
            );

      //double delta_phi_noise_driven(
      //      // -M_{phi} \frac{\delta F}{\delta \phi}
      //      const double& conc,
      //      const double& phi,
      //      const double& lap3dPhi,
      //      const double& TT,
      //      const std::vector<double>& hh
      //      );

      double dFdPhi(
            const double& conc,
            const double& Nv, // for conversion to/from walkers/conc.
            const double& phi,
            const double& TT,
            const std::vector<double>& hh
           );
            
      double c_p( const double& T); // heat capacity
      double latent_heat_c_p_inverse_integrand( 
                  const double& tt,// in case there's a time dependence
                  const double& TT,//domain variable of multiplicative func
                  const double& jump
                  );

}; // free_energy_gl class
} // SPF_NS

#endif
