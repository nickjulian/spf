
// File: free_energy_gl.cpp
//
#ifndef FREE_ENERGY_GL_CPP
#define FREE_ENERGY_GL_CPP
#include "free_energy_gl.hpp"

double SPF_NS::free_energy_gl::F(
   const double& conc, // in units of # walkers
   const double& Nv,  // to convert to concentration
   const double& phi,
   const double& TT,
   const double& lap1dPhiX,
   const double& lap1dPhiY,
   const double& lap1dPhiZ
   )
{
   double erfTemp;
   erfTemp = erf( ( TT- Fconstants[2]) / Fconstants[3]) +1;
   return (
         Fconstants[0] * phi * phi
         * pow( (conc/Nv) - (1 - Fconstants[1]*erfTemp), 2)
      ) + (
         Fconstants[0] * pow(1.0-phi,2)
         * pow( Fconstants[1]* erfTemp - (conc/Nv),2)
      ) + (
         Fconstants[4] * pow(1.0-phi,4)
      //) + (
      //   Fconstants[5] * (lap1dPhiX + lap1dPhiY + lap1dPhiZ)
      //   TODO: is the laplacian of phi inside F???
      );
}

double SPF_NS::free_energy_gl::dFdPhi(
            const double& conc, // in units of # of walkers
            const double& Nv, // for conversion to/from walkers/conc.
            const double& phi, // value of phi from previous timestep
            const double& TT,
            const std::vector<double>& hh
      )
{
   double erfTemp;
   erfTemp = erf( ( TT- Fconstants[2]) / Fconstants[3]) +1;

   return
      ( 2 * Fconstants[0]
         * phi
         * pow( (conc/Nv) - (1 - Fconstants[1] * erfTemp), 2)
      ) - (
         2 * Fconstants[0]
         * (1 - phi)
         * pow( Fconstants[1] * erfTemp - (conc/Nv), 2)
      ) - (
         4 * Fconstants[4] * pow( 1.0 - phi, 3)
      ); // NOTE: refer to Equation 4.16 page 33 of Nikolas Provatas and
         //  Ken Elder book 'Phase-Field Methods in Materials Science and
         //  Engineering'
         // There the laplacian of phi is inside an equation for
         //  \frac{\delta F}{\delta \phi} where \nabla^{2} \phi is added to
         //  \frac{\partial f}{\partial \phi} with an order energy W.
}

// double SPF_NS::free_energy_gl::delta_phi_noise_driven(
//            const double& conc,
//            const double& phi,
//            const double& lap3dPhi,
//            const double& TT,
//            const std::vector<double>& hh)
//{
//   TODO: add noise \sqrt{mu} \circ dB
//}

int SPF_NS::free_energy_gl::delta_phi_marcus(
            std::list<double>& deltaPhiSequence, // output
            const int& conc_local_flux_sum,
            const double& conc_initial,
            const double& Nv, // for conversion to/from walkers/conc.
            const double& phi,
            const double& lap3dPhi,
            const double& TT,
            const std::vector<double>& hh
      )
{
   // Equation 4.16 page 33 of Nikolas Provatas and Ken Elder book 
   //  'Phase-Field Methods in Materials Science and Engineering'

   deltaPhiSequence.clear();

   double phi_updated( phi);
   double avgMobility(
            (mobility_phi[0] +mobility_phi[1] +mobility_phi[2])/3.0
            );

   // Evaluate the sequence of deltaPhis cause by the sequence in concs.
   // Iterate from conc_initial to conc_initial + conc_local_flux_sum.
   if ( conc_local_flux_sum > 0)
   {
      for ( int ii=1; ii <= conc_local_flux_sum; ++ii)
      {
         deltaPhiSequence.push_back(
                  avgMobility
                   * (  orderEnergy * lap3dPhi
                        - dFdPhi(
                           (conc_initial + ii),
                           Nv,
                           phi_updated, // phi for conc of conc_inital+ii
                           TT,
                           hh
                           )
                     )
               );
         phi_updated += deltaPhiSequence.back(); 
      }
   }
   else if ( conc_local_flux_sum < 0)
   {
      for ( int ii=-1; ii >= conc_local_flux_sum; --ii)
      {
         deltaPhiSequence.push_back(
                  avgMobility
                   * (  orderEnergy * lap3dPhi
                        - dFdPhi(
                           (conc_initial + ii),
                           Nv,
                           phi_updated,
                           TT,
                           hh
                           )
                     )
               );
         phi_updated += deltaPhiSequence.back(); 
      }
   }
   return EXIT_SUCCESS;
}

double SPF_NS::free_energy_gl::delta_phi_noiseless(
            const double& conc, // in units of # of walkers
            const double& Nv, // for conversion to/from walkers/conc.
            const double& phi,
            const double& lap3dPhi,
            const double& TT,
            const std::vector<double>& hh
      )
{
   //std::cout << "dFdPhi( conc, Nv, phi, TT, hh): "
   //   << dFdPhi( conc, Nv, phi, TT, hh) << std::endl; // debug
   return (( mobility_phi[0] + mobility_phi[1] + mobility_phi[2])/3.0)
      * (  orderEnergy * lap3dPhi - dFdPhi( conc, Nv, phi, TT, hh));
   // Equation 4.16 page 33 of Nikolas Provatas and Ken Elder book 
   //  'Phase-Field Methods in Materials Science and Engineering'
}

double SPF_NS::free_energy_gl::delta_T_decoupled(
            //  \frac{\partial T}{\partial t} = D_{T} \nabla^{2} T
            const double& TT,
            const double& lapT, // laplacian of T
            const double& deltaPhi,
            const double& dt
            )
      {
         return thermalDiffusivity * lapT * dt;
      }

double SPF_NS::free_energy_gl::delta_T(
            //  \frac{\partial T}{\partial t} = D_{T} \nabla^{2} T
            const double& TT,
            const double& lapT, // laplacian of T
            const double& deltaPhi,
            const double& dt
            )
      {
         return thermalDiffusivity * lapT * dt
                  + (latentHeat/c_p(TT)) *deltaPhi;
      }

double SPF_NS::free_energy_gl::delta_T_marcus(
            //  \frac{\partial T}{\partial t} = D_{T} \nabla^{2} T
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
            )
      {
         // TODO: iterate over deltaPhis and implement the marcus integral

         double driftDeltaT(thermalDiffusivity * lapT * dt); // to be added
         double jump_T( TT);

         // for each jump in phi, calculate the jump in TT using the marcus
         //  method 
         for ( std::list<double>::const_iterator
                  itrPhi=deltaPhiSequence.begin();
               itrPhi != deltaPhiSequence.end();
               ++itrPhi)
         {
            // integrate over the unit interval using the marcus method
            for( std::vector<double>::const_iterator
                     itrI=unitInterval.begin();
                  itrI != unitInterval.end();
                  ++itrI)
            {
               jump_T // jump destination
                  = marcusRK4(
                        //latent_heat_c_p_inverse,
                        //latent_heat_c_p_inverse_integrand,
                                    rk_dt,
                                    jump_T, // evolves with integration
                                    *itrI, // t0
                                    *itrPhi // jump in driving noise
                                    );
            }
            // Now jump_T is the value of TT after one of many jumps in Phi
            // Although the formulae contain more terms, I suspect they
            //  cancel in this case and are there only to maintain notation
            //  which represents the limit as dt goes to 0.
         }

         //return thermalDiffusivity * lapT * dt
         //         + (latentHeat/c_p(TT)) *deltaPhi;

         // integrate latentHeat/c_p(TT)
         // latent_heat_c_p_inverse( TT) // integrand

         //marcusIntegral(
         //   //double (*marcusIntegrand)(const double& tt,
         //   //                     const double& yy,
         //   //                     const double& dP), 
         //   //const double& rk_dt,    // increment of unit interval
         //   //const double& y0,    // path.back()
         //   //const double& t0,    // t_n, t_{n+1}=t_n + dt
         //   //const double& dP,// Y_{k} size of noise jump
         //   //double& jumpDestination
         //   );

         return driftDeltaT + (jump_T - TT);
      }

double SPF_NS::free_energy_gl::marcusRK4(
            const double& dt,
            const double& y0,
            const double& t0,
            const double& jump
            )
{
   // This function was brought into this class because using
   //  SPF_NS::marcusRK4() produced the error:
   //  "invalid use of non-static member function," referring to
   //  passing latent_heat_c_p_inverse_integrand to it.
   double k1, k2, k3, k4; 
   k1 = dt* latent_heat_c_p_inverse_integrand( t0, y0, jump); // == y0*jump
   k2 = dt* latent_heat_c_p_inverse_integrand( t0 + 0.5*dt, y0 + 0.5*k1, jump);
   k3 = dt* latent_heat_c_p_inverse_integrand( t0 + 0.5*dt, y0 + 0.5*k2, jump);
   k4 = dt* latent_heat_c_p_inverse_integrand( t0 + dt, y0 + k3, jump);

   return y0 + ONESIXTH *(k1 + 2.0*k2 + 2.0*k3 + k4 );
}

double SPF_NS::free_energy_gl::latent_heat_c_p_inverse_integrand(
      const double& tt, // in case there's a time dependence (there isn't)
      const double& TT, // domain variable of multiplicative function
      const double& jump // Y_k // jump in driving noise
      )
   {
      return jump * latentHeat / c_p( TT);
   }

double SPF_NS::free_energy_gl::c_p( const double& TT)
{
   // Heat capacity
   // c_{p}(T,\alpha,\beta,\gamma,cbase) = -0.05 \frac{\alpha^{2}}{T^{2} K_{eq}(T,\alpha,\gamma)} (0.5-\frac{2+ \frac{1}{K_{eq}(T,\alpha,\gamma)}}{2 \sqrt{(2+\frac{1}{K_{eq}(T,\alpha,\gamma)})^{2}-4}} + cbase
   // K_{eq}(T,\alpha,T0) = \exp(\frac{-1*alpha}{T} - T0)
   double K_eq;
   K_eq = exp( (-1*alpha/TT) +T0);
   return 
      prefactor
      * ((alpha*alpha)/(TT*TT * K_eq))
      * (0.5 - (2+(1/K_eq))/(2 * sqrt(pow(2+(1/K_eq),2)-4)))
      + cbase;
}

double SPF_NS::free_energy_gl::dF_dconc( 
            const double& conc, // in units of # of walkers
            const double& Nv, // for conversion to/from walkers/conc.
            const double& phi,
            const double& TT
            )
{
   // dF/dconc
   // = 2 c1 (1-phi)^2 (conc-(1-c2(erf((T-c3)/c4)+1)))
   //   - 2 c1 phi^2 (c2(1+erf((T-c3)/c4)) - conc)
   double erfTemp;
   erfTemp = erf( ( TT- Fconstants[2]) / Fconstants[3]) +1; // debug
   return 
      fmax(2 * Fconstants[0] * pow(1 - phi, 2)
       * ((conc/Nv) - (1 - Fconstants[1] * erfTemp))
      - 2 * Fconstants[0] * phi * phi
       * (Fconstants[1] * erfTemp - (conc/Nv)), 0.0);// output positive calc or 0
}

double SPF_NS::free_energy_gl::entropy( 
            const double& conc, // in units of # of walkers
            const double& phi
            )
{
   return -1.0 * kb * ( conc * log( conc) +  (1.0-phi) * log( 1.0 - phi));
}

#endif
