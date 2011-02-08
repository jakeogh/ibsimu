/*! \file types.hpp
 *  \brief Base types
 */

/* Copyright (c) 2005-2010 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * taneli.kalvas@jyu.fi.
 * 
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
 */

#ifndef TYPES_HPP
#define TYPES_HPP 1


/*! \brief Geometry mode enum.
 *
 *  Simulation geometry can be either 1D (\a MODE_1D), planar 2D (\a
 *  MODE_2D), planar 3D (\a MODE_3D) or it can be defined as cylindrical
 *  symmetrical 2D with coordinates \a x and \a r (\a MODE_CYL).
 *
 *  The geometry mode selects the active axes for calculation and
 *  fields. For MODE_1D only x-axis (axis 0) is active. For MODE_2D
 *  the x- and y-axes are active (axes 0 and 1). For MODE_CYL
 *  the x- and r-axes are active (axes 0 and 1). For MODE_3D
 *  the x- and y- and z-axes are active (axes 0, 1 and 2).
 */
enum geom_mode_e {
    MODE_1D = 0,
    MODE_2D,
    MODE_CYL, 
    MODE_3D
};


/*! \brief Field extrapolation mode
 *
 *  This parameter is used to control the behaviour of the field
 *  evaluators outside the defined area. The field value can be
 *  extrapolated from the closest defined points (\a
 *  FIELD_EXTRAPOLATE), the field can be mirrored (\a FIELD_MIRROR),
 *  the field evaluator can simply return zero (\a FIELD_ZERO) or the
 *  field evaluator can return not-a-number, NaN (\a FIELD_NAN ).
 */
enum field_extrpl_e {
    FIELD_EXTRAPOLATE = 0,
    FIELD_MIRROR,
    FIELD_ZERO,
    FIELD_NAN
};


/*! \brief Field type.
 *
 *  Indicator for field type.
 */
enum field_type_e {
    FIELD_NONE = 0,
    FIELD_EPOT,
    FIELD_SCHARGE,
    FIELD_TRAJDENS,
    FIELD_EFIELD,
    FIELD_EFIELD_X,
    FIELD_EFIELD_Y,
    FIELD_EFIELD_Z,
    FIELD_BFIELD,
    FIELD_BFIELD_X,
    FIELD_BFIELD_Y,
    FIELD_BFIELD_Z
};


/*! \brief Field diagnostic type.
 *
 *  \deprecated Provided for compatibility only. Replaced by
 *  multipurpose field type selector field_type_e.
 */
typedef field_type_e field_diag_type_e;


/*! \brief Boundary type.
 *
 *  Boundary conditions for solids and simulation box boundaries. See
 *  class Bound for more information.
 */
enum bound_e {
    BOUND_DIRICHLET = 0, /*!< \brief Dirichlet boundary condition */
    BOUND_NEUMANN        /*!< \brief Neumann (or natural) boundary condition */
};


/*! \brief Coordinate axis identifier.
 */
enum coordinate_axis_e {
    AXIS_X = 0, /*!< \brief X axis */
    AXIS_Y,     /*!< \brief Y axis */
    AXIS_R,     /*!< \brief R axis */
    AXIS_Z      /*!< \brief Z axis */
};


/*! \brief String describing axis names without unit.
 *
 *  Contains strings: "x", "y", "r" and "z".
 */
extern const char *coordinate_axis_string[];


/*! \brief String describing axis names with unit.
 *
 *  Contains strings: "x (m)", "y (m)", "r (m)" and "z (m)".
 */
extern const char *coordinate_axis_string_with_unit[];


/*! \brief Type of diagnostic for trajectories.
 *
 *  O-, P- and Q-axes are diagnostic axes defined by user.
 */
enum trajectory_diagnostic_e {
    DIAG_NONE = 0, /*!< \brief Dummy diagnostic. Does nothing. */
    DIAG_T,        /*!< \brief Time (s) */
    DIAG_X,        /*!< \brief X-axis position (m) */
    DIAG_VX,       /*!< \brief X-axis velocity (m/s) */
    DIAG_Y,        /*!< \brief Y-axis position (m) */
    DIAG_R,        /*!< \brief Radial position (m) */
    DIAG_VY,       /*!< \brief Y-axis velocity (m/s) */
    DIAG_VR,       /*!< \brief Radial velocity (m/s) */
    DIAG_W,        /*!< \brief Angular velocity (rad/s) */
    DIAG_VTHETA,   /*!< \brief Tangential velocity (m/s) */
    DIAG_Z,        /*!< \brief Z-axis position (m) */
    DIAG_VZ,       /*!< \brief Z-axis velocity (m/s) */
    DIAG_O,        /*!< \brief O-axis position (m) */
    DIAG_VO,       /*!< \brief O-axis velocity (m/s) */
    DIAG_P,        /*!< \brief P-axis position (m) */
    DIAG_VP,       /*!< \brief P-axis velocity (m/s) */
    DIAG_Q,        /*!< \brief Q-axis position (m) */
    DIAG_VQ,       /*!< \brief Q-axis velocity (m/s) */
    DIAG_XP,       /*!< \brief \f$v_x/v_q\f$, where direction q is normal to diagnostic plane (rad) */
    DIAG_YP,       /*!< \brief \f$v_y/v_q\f$, where direction q is normal to diagnostic plane (rad) */
    DIAG_RP,       /*!< \brief \f$v_r/v_q\f$, where direction q is normal to diagnostic plane (rad) */
    DIAG_AP,       /*!< \brief \f$v_{\theta}/v_q\f$, where direction q is normal to diagnostic plane (rad) */
    DIAG_ZP,       /*!< \brief \f$v_z/v_q\f$, where direction q is normal to diagnostic plane (rad) */
    DIAG_OP,       /*!< \brief \f$v_o/v_q\f$, where direction q is normal to diagnostic plane (rad) */
    DIAG_PP,       /*!< \brief \f$v_p/v_q\f$, where direction q is normal to diagnostic plane (rad) */
    DIAG_CURR,     /*!< \brief Current (I) */
    DIAG_EK,       /*!< \brief Kinetic energy (J) */
    DIAG_QM        /*!< \brief Charge per mass (C/g) */
};


/*! \brief String describing diagnostic without unit.
 *
 *  Contains strings: "none", "t", "x", "v_x", "y", ... Greek letters
 *  are typed with LaTeX notation for correct output in plots.
 */
extern const char *trajectory_diagnostic_string[];


/*! \brief String describing diagnostic with unit.
 *
 *  Contains strings: "none ()", "t (s)", "x (m)", "v_x (m/2)", "y
 *  (m)", ... Greek letters are typed with LaTeX notation for correct
 *  output in plots.
 */
extern const char *trajectory_diagnostic_string_with_unit[];


#endif
