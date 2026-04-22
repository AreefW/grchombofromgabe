/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COSMOMATTERCCZ4RHS_HPP_
#define COSMOMATTERCCZ4RHS_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4RHS.hpp"
#include "Cell.hpp"
#include "CosmoMovingPunctureGauge.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//!  Calculates RHS using CCZ4 including matter terms, and matter variable
//!  evolution, with runtime cosmological gauge K_mean support.
/*!
     This class mirrors MatterCCZ4RHS but fixes the ScalarFieldCosmo use-case
   where K_mean must be updated at runtime for CosmoMovingPunctureGauge.
*/

template <class matter_t, class deriv_t = FourthOrderDerivatives>
class CosmoMatterCCZ4RHS : public CCZ4RHS<CosmoMovingPunctureGauge, deriv_t>
{
  public:
    using gauge_t = CosmoMovingPunctureGauge;
    using CCZ4 = CCZ4RHS<gauge_t, deriv_t>;

    using params_t = CCZ4_params_t<typename gauge_t::params_t>;

    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    template <class data_t>
    using MatterDiff2Vars = typename matter_t::template Diff2Vars<data_t>;

    template <class data_t>
    using CCZ4Vars = typename CCZ4::template Vars<data_t>;

    template <class data_t>
    using CCZ4Diff2Vars = typename CCZ4::template Diff2Vars<data_t>;

    // Inherit the variable definitions from CCZ4RHS + matter_t
    template <class data_t>
    struct Vars : public CCZ4Vars<data_t>, public MatterVars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            CCZ4Vars<data_t>::enum_mapping(mapping_function);
            MatterVars<data_t>::enum_mapping(mapping_function);
        }
    };

    template <class data_t>
    struct Diff2Vars : public CCZ4Diff2Vars<data_t>,
                       public MatterDiff2Vars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            CCZ4Diff2Vars<data_t>::enum_mapping(mapping_function);
            MatterDiff2Vars<data_t>::enum_mapping(mapping_function);
        }
    };

    //! Constructor
    CosmoMatterCCZ4RHS(matter_t a_matter, params_t a_params, double a_dx,
                       double a_sigma, int a_formulation = CCZ4RHS<>::USE_CCZ4,
                       double a_G_Newton = 1.0, double a_K_mean = 0.0);

    //! Update K_mean at runtime (call per RHS evaluation if needed)
    void set_K_mean(double a_K_mean);

    //! The compute member which calculates the RHS at each point in the box
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    //! Add EM Tensor terms to CCZ4 rhs
    template <class data_t>
    void add_emtensor_rhs(
        Vars<data_t> &matter_rhs,         //!< RHS data at this point
        const Vars<data_t> &vars,         //!< variable values at this point
        const Vars<Tensor<1, data_t>> &d1 //!< first derivatives at this point
    ) const;

    // Class members
    matter_t my_matter;      //!< The matter object, e.g. a scalar field.
    const double m_G_Newton; //!< Newton's constant, set to one by default.
    double m_K_mean;         //!< Runtime K_mean used by cosmological gauge.
};

#include "CosmoMatterCCZ4RHS.impl.hpp"

#endif /* COSMOMATTERCCZ4RHS_HPP_ */