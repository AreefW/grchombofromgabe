/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDGRIDSTAGGINGCRITERION_HPP_
#define FIXEDGRIDSTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

class FixedGridsTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const int m_level;
    const int m_max_level;
    const std::array<double, CH_SPACEDIM> m_center;

  public:
    FixedGridsTaggingCriterion(const double dx, const int a_level,
                               const int a_max_level, const double a_L,
                               const std::array<double, CH_SPACEDIM> a_center)
        : m_dx(dx), m_L(a_L), m_level(a_level), m_max_level(a_max_level),
          m_center(a_center){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        data_t criterion = 0.0;
        // make sure the inner part is regridded around the horizon
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4
        // double ratio = pow(2.0, -(m_level + 2.0));
        // int level_diff = m_max_level - m_level;
        // // pout() << "Lelev_diff = " << level_diff << endl;
        // double ratio = 0.0;
        // if (level_diff <= m_max_level){
        //   ratio = (1. + 0.1 * level_diff);
        // }
        // else {
        //   ratio = 0.;
        // }

        double ratio = 0.0;

        if (m_level < m_max_level)
        {
            ratio = 1.0 + 0.125 * (m_max_level - m_level);
        }
        else
        {
            ratio = 1.0;
        }
        // if (m_level == 0){
        //   ratio = 1.4;
        // }
        // else if (m_level == 1){
        //   ratio = 1.2;
        // }
        // else if (m_level == 2){
        //   ratio = 1.;
        // }
        // else {
        //   ratio = 0.;
        // }
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        // const data_t max_abs_xy = simd_max(abs(coords.x), abs(coords.y));
        // const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));
        // auto in_box_x = simd_compare_lt(abs(coords.x - m_center[0]), m_L);
        // auto abs_xy = simd_compare_lt(in_box_x, abs(coords.x - m_center[0]),
        // m_L); auto in_box_y = simd_compare_lt(abs(coords.y) - abs(coords.x),
        // m_L); auto x_distance = simd_conditional(inbox_x, abs(coords.x -
        // m_center[0]), ); auto inbox_x = simd_compare_lt(coords.x -
        // m_center[0], m_L); auto inbox_y = simd_compare_lt(coords.y -
        // m_center[1], m_L); auto inbox_z = simd_compare_lt(coords.z -
        // m_center[2], m_L); auto inbox = simd_compare_lt(abs(coords.x -
        // m_center[0]) + abs(coords.y - m_center[1]) + abs(coords.z -
        // m_center[2]), (3 * m_L * ratio)); auto regrid = inbox_x * inbox_y *
        // inbox_z;

        // const data_t r = coords.get_radius();
        // auto regrid = simd_compare_lt(r, m_L * ratio);
        // pout() << "abs_x = " << coords.x << "-" << m_center[0] << " = " <<
        // abs(coords.x - m_center[0]) << endl;

        // const data_t max_abs_xy = simd_max(abs(coords.x - m_center[0]),
        // abs(coords.y - m_center[1])); const data_t max_abs_xyz =
        // simd_max(max_abs_xy, abs(coords.z - m_center[2])); auto regrid =
        // simd_compare_lt(max_abs_xyz, m_L * ratio);

        const data_t max_abs_xy = simd_max(abs(coords.x), abs(coords.y));
        const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));
        auto regrid = simd_compare_lt(max_abs_xyz, m_L * ratio);
        criterion = simd_conditional(regrid, 100.0, criterion);

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* FIXEDGRIDSTAGGINGCRITERION_HPP_ */
