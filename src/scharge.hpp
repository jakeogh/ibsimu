/*! \file scharge.hpp
 *  \brief Header file defining space charge.
 */


#ifndef SCHARGE_HPP
#define SCHARGE_HPP 1


#include "scalarfield.hpp"
#include "particles.hpp"



/*! \page p_scharge Space charge
 *
 *  Space charge in simulations is handled as space charge \b density
 *  (C/m3) stored in a ScalarField object.
 *
 *  The space charge field is cleared by the particle iterator before
 *  calculating particle trajectories. While calculating trajectories
 *  the iterator calls scharge_add_from_trajectory() at each mesh
 *  location the trajectory passes through. The charge deposited by
 *  the trajetory to the mesh is calculated. At this phase the \a
 *  scharge field contains charge per unit length (C/m) in case of 2D.
 *  and charge (C) in other cases. The space charge map is converted
 *  to space charge \b density map by a call to scharge_finalize() by
 *  the particle iterator at the end of particle calculation.
 */


/*! \brief Finalize space charge calculation.
 *
 *  Converts charge map built by calls to
 *  scharge_add_from_trajectory() to space charge density map suitable
 *  for potential solvers. Also corrects boundary space charge values.
 *
 */
void scharge_finalize( ScalarField &scharge );


/*! \brief Function for adding charge to space charge density map from
 *  particle trajectory in 2d simulation.
 *
 *  Adds contribution to space charge density map \a scharge from a
 *  particle carrying current per length \a IQ travelling from \a x1
 *  to \a x2. Charge density is added to the mesh assuming that it is
 *  localized at the midpoint of \a x1 and \a x2. This function is
 *  supposed to be used so that \a x1 and \a x2 are at the mesh
 *  intersection points, through which particle trajectory has passed.
 */
void scharge_add_from_trajectory( ScalarField &scharge, double IQ, 
				  const ParticleP2D &x1, const ParticleP2D &x2 );


/*! \brief Function for adding charge to space charge density map from
 *  particle trajectory in cylindrically symmetric simulation.
 *
 *  Same as for scharge_add_from_trajectory(). Now IQ is real current
 *  (A).
 */
void scharge_add_from_trajectory( ScalarField &scharge, double IQ, 
				  const ParticlePCyl &x1, const ParticlePCyl &x2 );


/*! \brief Function for adding charge to space charge density map from
 *  particle trajectory in 3d simulation.
 *
 *  Same as for scharge_add_from_trajectory(). Now IQ is real current
 *  (A).
 */
void scharge_add_from_trajectory( ScalarField &scharge, double IQ, 
				  const ParticleP3D &x1, const ParticleP3D &x2 );


#endif













