/*! \file epot_efield.hpp
 *  \brief Header file defining electric field class based on
 *  on-line interpolation of electric potential.
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
 *  \a x (9 points in 2D and 27 points in 3D).  A second order
 *  polynomial is fitted to these points to get a linear (continuous)
 *  interpolation of electric field. Use of solid mesh points is
 *  avoided in the interpolation because they would mess with the
 *  virtual edge points.
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













