/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SCALARFIELD_HPP_)
#error "This file should only be included through ScalarField.hpp"
#endif

#ifndef SCALARFIELD_IMPL_HPP_
#define SCALARFIELD_IMPL_HPP_

// Two-field kinetic preheating system (arXiv:2511.02059)
// Lagrangian: L = -1/2(d phi)^2 - W(phi)/2 (d chi)^2 - V(phi)
// where W(phi) = exp(2 phi / Mu) is the kinetic coupling and
// chi is the axion field (called 'axion' in code to avoid clash with
// GRChombo's conformal factor 'chi').
//
// Canonical momenta (Eqs. 4):
//   Pi    = (1/lapse)(d_t phi  - beta^i d_i phi)
//   Theta = (1/lapse)(d_t axion - beta^i d_i axion)

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> ScalarField<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;

    // call the function which computes the em tensor excluding the potential
    emtensor_excl_potential(out, vars, d1, h_UU, chris_ULL);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi = 0.0;

    // compute potential and add contributions to EM Tensor
    my_potential.compute_potential(V_of_phi, dVdphi, vars);

    out.rho += V_of_phi;
    out.S += -3.0 * V_of_phi;
    FOR(i, j) { out.Sij[i][j] += -vars.h[i][j] * V_of_phi / vars.chi; }

    return out;
}

// Calculate the stress energy tensor elements excluding the potential
// For Lagrangian L = -1/2(d phi)^2 - W(phi)/2 (d axion)^2 - V(phi):
//   T_mu_nu = d_mu phi d_nu phi + W d_mu axion d_nu axion
//             - g_mu_nu (1/2(d phi)^2 + W/2(d axion)^2 + V)
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void ScalarField<potential_t>::emtensor_excl_potential(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &h_UU,
    const Tensor<3, data_t> &chris_ULL) const
{
    // Get kinetic coupling W(phi) = exp(2 phi / Mu)
    data_t W_of_phi = 0.0;
    data_t dWdphi = 0.0;
    my_potential.compute_kinetic_coupling(W_of_phi, dWdphi, vars);

    // Vt_phi = g^{mu nu} d_mu phi d_nu phi = -Pi^2 + chi * h^{ij} d_i phi d_j phi
    data_t Vt_phi = -vars.Pi * vars.Pi;
    FOR(i, j) { Vt_phi += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j]; }

    // Vt_axion = g^{mu nu} d_mu axion d_nu axion (without W factor)
    data_t Vt_axion = -vars.Pi_axion * vars.Pi_axion;
    FOR(i, j)
    {
        Vt_axion += vars.chi * h_UU[i][j] * d1.axion[i] * d1.axion[j];
    }

    // Combined kinetic term: (d phi)^2 + W (d axion)^2
    data_t Vt_total = Vt_phi + W_of_phi * Vt_axion;

    // Calculate components of EM Tensor
    // S_ij = T_ij (spatial stress)
    FOR(i, j)
    {
        out.Sij[i][j] = -0.5 * vars.h[i][j] * Vt_total / vars.chi +
                        d1.phi[i] * d1.phi[j] +
                        W_of_phi * d1.axion[i] * d1.axion[j];
    }

    // S = Tr S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    // S_i (lower index) = -n^a T_{ai}
    FOR(i)
    {
        out.Si[i] =
            -d1.phi[i] * vars.Pi - W_of_phi * d1.axion[i] * vars.Pi_axion;
    }

    // rho = n^a n^b T_{ab}
    out.rho = vars.Pi * vars.Pi + W_of_phi * vars.Pi_axion * vars.Pi_axion +
              0.5 * Vt_total;
}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void ScalarField<potential_t>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // call the function for the rhs excluding the potential
    matter_rhs_excl_potential(total_rhs, vars, d1, d2, advec);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi = 0.0;
    my_potential.compute_potential(V_of_phi, dVdphi, vars);

    // adjust RHS for the potential term in the inflaton equation (Eq. 5)
    // d_t Pi contribution: -alpha * dV/dphi
    total_rhs.Pi += -vars.lapse * dVdphi;
    // axion has no potential coupling (V = V(phi) only)
}

// The RHS excluding the potential terms
// Implements Eqs. (5) and (6) of arXiv:2511.02059, extended with BSSN
// conformal decomposition used in GRChombo (chi = e^{-4 phi_conf}).
//
// Inflaton Pi equation (Eq. 5):
//   d_0 Pi = beta^k d_k Pi + alpha K Pi
//            + gamma^{ij}(d_i alpha d_j phi + alpha d_i d_j phi
//                         - alpha Gamma^k_{ij} d_k phi)
//            - alpha/2 * dW/dphi * (gamma^{ij} d_i axion d_j axion - Theta^2)
//   [potential term -alpha dV/dphi added in add_matter_rhs]
//
// Axion Theta equation (Eq. 6):
//   d_0 Theta = beta^k d_k Theta + alpha K Theta
//               + gamma^{ij}(d_i alpha d_j axion + alpha d_i d_j axion
//                            - alpha Gamma^k_{ij} d_k axion)
//               + alpha/W * dW/dphi * (gamma^{ij} d_i phi d_j axion - Pi*Theta)
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void ScalarField<potential_t>::matter_rhs_excl_potential(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // Get kinetic coupling W(phi) and dW/dphi
    data_t W_of_phi = 0.0;
    data_t dWdphi = 0.0;
    my_potential.compute_kinetic_coupling(W_of_phi, dWdphi, vars);

    // Vt_axion = g^{mu nu} d_mu axion d_nu axion (without W)
    // = -Theta^2 + chi * h^{ij} d_i axion d_j axion
    data_t Vt_axion = -vars.Pi_axion * vars.Pi_axion;
    FOR(i, j)
    {
        Vt_axion += vars.chi * h_UU[i][j] * d1.axion[i] * d1.axion[j];
    }

    // -------------------------------------------------------------------------
    // Inflaton phi evolution: d_t phi = alpha Pi + beta^i d_i phi
    // -------------------------------------------------------------------------
    rhs.phi = vars.lapse * vars.Pi + advec.phi;

    // -------------------------------------------------------------------------
    // Inflaton Pi evolution (Eq. 5, standard KG terms)
    // -------------------------------------------------------------------------
    rhs.Pi = vars.lapse * vars.K * vars.Pi + advec.Pi;

    FOR(i, j)
    {
        // Standard KG terms in conformal BSSN (conformal factor chi):
        //   chi * h^{ij} (alpha d_i d_j phi + d_i alpha d_j phi)
        //   - alpha/2 * h^{ij} d_j chi d_i phi   [conformal correction to Gamma]
        rhs.Pi += h_UU[i][j] * (-0.5 * d1.chi[j] * vars.lapse * d1.phi[i] +
                                vars.chi * vars.lapse * d2.phi[i][j] +
                                vars.chi * d1.lapse[i] * d1.phi[j]);
        FOR(k)
        {
            // -alpha chi h^{ij} Gamma^k_{ij}[h] d_k phi
            rhs.Pi += -vars.chi * vars.lapse * h_UU[i][j] *
                      chris.ULL[k][i][j] * d1.phi[k];
        }
    }

    // Kinetic coupling term for Pi (Eq. 5):
    //   -alpha/2 * dW/dphi * (gamma^{ij} d_i axion d_j axion - Theta^2)
    //   = -alpha/2 * dW/dphi * Vt_axion
    rhs.Pi += -vars.lapse * 0.5 * dWdphi * Vt_axion;

    // -------------------------------------------------------------------------
    // Axion field evolution: d_t axion = alpha Theta + beta^i d_i axion
    // -------------------------------------------------------------------------
    rhs.axion = vars.lapse * vars.Pi_axion + advec.axion;

    // -------------------------------------------------------------------------
    // Axion Theta evolution (Eq. 6, standard KG terms for the axion)
    // -------------------------------------------------------------------------
    rhs.Pi_axion = vars.lapse * vars.K * vars.Pi_axion + advec.Pi_axion;

    FOR(i, j)
    {
        // Standard KG terms in conformal BSSN applied to the axion field:
        rhs.Pi_axion +=
            h_UU[i][j] * (-0.5 * d1.chi[j] * vars.lapse * d1.axion[i] +
                          vars.chi * vars.lapse * d2.axion[i][j] +
                          vars.chi * d1.lapse[i] * d1.axion[j]);
        FOR(k)
        {
            rhs.Pi_axion += -vars.chi * vars.lapse * h_UU[i][j] *
                         chris.ULL[k][i][j] * d1.axion[k];
        }
    }

    // Kinetic coupling term for Theta (Eq. 6):
    //   + alpha/W * dW/dphi * (gamma^{ij} d_i phi d_j axion - Pi * Theta)
    //   = alpha * (dW/dphi / W) * (chi h^{ij} d_i phi d_j axion - Pi * Theta)
    data_t phi_axion_inner = -vars.Pi * vars.Pi_axion;
    FOR(i, j)
    {
        phi_axion_inner +=
            vars.chi * h_UU[i][j] * d1.phi[i] * d1.axion[j];
    }
    rhs.Pi_axion += vars.lapse * (dWdphi / W_of_phi) * phi_axion_inner;
}

#endif /* SCALARFIELD_IMPL_HPP_ */
