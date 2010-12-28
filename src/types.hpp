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
 *  FIELD_EXTRAPOLATE), the field can be mirrored (\a FIELD_MIRROR) or
 *  the field evaluator can simply return zero (\a FIELD_ZERO).
 */
enum field_extrpl_e {
    FIELD_EXTRAPOLATE = 0,
    FIELD_MIRROR,
    FIELD_ZERO
};


/*! \brief Boundary type enum.
 *
 *  %Solid boundary may have either essential, Dirichlet boundary
 *  condition \a BOUND_DIRICHLET or natural, Neumann boundary
 *  condition \a BOUND_NEUMANN.
 */
enum bound_e {
    BOUND_DIRICHLET = 0,
    BOUND_NEUMANN
};


/*! \brief Coordinate axis identifier.
 */
enum coordinate_axis_e {
    AXIS_X = 0,
    AXIS_Y,
    AXIS_R,
    AXIS_Z
};


extern const char *coordinate_axis_string[];
extern const char *coordinate_axis_string_with_unit[];


/*! \brief Type of diagnostic for trajectories.
 */
enum trajectory_diagnostic_e {
    DIAG_NONE = 0,
    DIAG_T,
    DIAG_X,
    DIAG_VX,
    DIAG_Y,
    DIAG_R,
    DIAG_VY,
    DIAG_VR,
    DIAG_W,
    DIAG_VTHETA,
    DIAG_Z,
    DIAG_VZ,
    DIAG_CURR,
    DIAG_XP,
    DIAG_YP,
    DIAG_RP,
    DIAG_AP,
    DIAG_ZP,
    DIAG_EK,
    DIAG_QM
};


extern const char *trajectory_diagnostic_string[];
extern const char *trajectory_diagnostic_string_with_unit[];






#endif


















