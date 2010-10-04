/*! \file epot_efield.hpp
 *  \brief Header file for epot_efield.hpp
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
 * tvkalvas@cc.jyu.fi.
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

#ifndef EPOT_EFIELD_HPP
#define EPOT_EFIELD_HPP 1


#include "efield.hpp"
#include "geometry.hpp"
#include "scalarfield.hpp"
#include "vec3d.hpp"


enum efield_extrpl_e {
    EFIELD_EXTRAPOLATE = 0,
    EFIELD_MIRROR,
    EFIELD_ZERO
};


/*! \brief Electric field class implementation based on on-line
 *  interpolation of electric potential.
 *
 *  %EpotEfield contains pointers to Geometry and to ScalarField
 *  electric potential (epot). In 1D %EpotEfield uses three closest
 *  mesh neighbours of point \a x to calculate the electric field at
 *  \a x in 1D (9 points in 2D and 27 points in 3D).  A second order
 *  polynomial is fitted to the potential field to get a linear
 *  (continuous) interpolation of electric field. Use of solid mesh
 *  points is avoided in the interpolation because they would mess
 *  with the virtual edge points.
 *
 *  The behaviour of the interpolation function outside mesh points
 *  can be programmed with set_extrapolation() function. Behaviour
 *  defaults to extrapolation using closest electric potential points. 
 */
class EpotEfield : public Efield {

    efield_extrpl_e      _extrpl[6];   /*!< \brief What to return outside geometry. */

    const Geometry      &_g;           /*!< \brief Reference to geometry. */
    const ScalarField   &_epot;        /*!< \brief Reference to electric potential. */

public:

    /*! \brief Constructor.
     */
    EpotEfield( const Geometry &g, const ScalarField &epot );

    /*! \brief Destructor.
     */
    ~EpotEfield() {}

    /*! \brief Set the behaviour of electric field interpolation
     *  outside mesh points (extrapolation).
     *
     *  The interpolation function behaviour can be set separately for
     *  each boundary. This is done by setting the desired properties
     *  to the \a extrpl array. The interpolation function can use an
     *  extrapolation of the last three electric potential values for
     *  calculation of electric field (EFIELD_EXTRAPOLATE) or it can
     *  return the mirror of the electric field across the mesh
     *  boundary (EFIELD_MIRROR) or it can return a zero electric
     *  field outside the mesh.
     *
     *  The use of EFIELD_MIRROR in case of symmetric cases, where
     *  beam is traversing next to the geometry boundary, is necessary
     *  to get physical results.
     *
     *  Very far (double the size of the simulation box) the field
     *  evaluator will always return zero.
     */
    void set_extrapolation( efield_extrpl_e extrpl[6] ) {
	memcpy( _extrpl, extrpl, 6*sizeof(efield_extrpl_e) );
    }

    /*! \brief Operator for getting interpolated electric field value
     *  at \a x.
     */
    Vec3D operator()( Vec3D x ) const;
};


#endif















