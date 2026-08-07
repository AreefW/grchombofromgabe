/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"

class Potential
{
  public:
    struct params_t
    {
        double scalar_mass;
    };

  private:
    params_t m_params;
    static constexpr double Mu = 0.0468; // The value Eve calculated IF CHANGING MUST RECOMPILE!!!!

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi
       // V_of_phi = 0.5 * pow(m_params.scalar_mass * vars.phi, 2.0); //m^2 * mu^2 /2 (1-exp(-phi/mu))^2
	V_of_phi = 0.5 * pow(m_params.scalar_mass*Mu*(1-exp(-vars.phi/Mu)), 2.0);
        // The potential gradient at phi
       // dVdphi = pow(m_params.scalar_mass, 2.0) * vars.phi;
	dVdphi = pow(m_params.scalar_mass, 2.0) * Mu*exp(-vars.phi/Mu)*(1-exp(-vars.phi/Mu));
    }

    //! Set the kinetic coupling function W(phi) = exp(2*phi/Mu)
    //! This encodes the kinetic preheating coupling to the axion field chi
    //! (arXiv:2511.02059, Eq. 1: L = -1/2(d phi)^2 - W(phi)/2 (d chi)^2 - V(phi))
    template <class data_t, template <typename> class vars_t>
    void compute_kinetic_coupling(data_t &W_of_phi, data_t &dWdphi,
                                  const vars_t<data_t> &vars) const
    {
        // W(phi) = exp(2 phi / Mu)
        W_of_phi = exp(2.0 * vars.phi / Mu);
        // dW/dphi = (2/Mu) * W
        dWdphi = (2.0 / Mu) * W_of_phi;
    }
};

#endif /* POTENTIAL_HPP_ */

