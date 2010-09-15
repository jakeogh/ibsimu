/*! \file particles.cpp
 *  \brief Source code for particles.cpp
 */

/* Copyright (c) 2005-2009 Taneli Kalvas. All rights reserved.
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

#include "particles.hpp"
#include "trajectory.hpp"
#include <iostream>
#include <iomanip>


//#define DEBUG_PARTICLE_DERIVATIVES 1


int ParticleP2D::get_derivatives( double t, const double *x, double *dxdt, void *data )
{
    Vec3D E, B, xc( x[0], x[2], 0.0 );
    ParticleIteratorData *pi = (ParticleIteratorData *)data;
    
    if( pi->_efield )
	E = (*pi->_efield)( xc );
    if( pi->_bfield )
	B = (*pi->_bfield)( xc );
    //std::cout << "E=(" << E << ") at x=(" << xc << ")\n";
    
    /* Positions: dx/dt = vx, dy/dt = vy */
    dxdt[0] = x[1];
    dxdt[2] = x[3];
    
    /* Velocities dvx/dt = ax, dvy/dt = ay */
    dxdt[1] = pi->_qm * (E[0] + x[3]*B[2]);
    dxdt[3] = pi->_qm * (E[1] - x[1]*B[2]);

    //std::cout << "dxdt=(" 
    //<< dxdt[0] << " "
    //<< dxdt[1] << " "
    //<< dxdt[2] << " "
    //<< dxdt[3] << ")\n";
	
    return( GSL_SUCCESS );
}


int ParticleP2D::trajectory_intersections_at_plane( std::vector<ParticleP2D> &intsc, 
						    int crd, double val,
						    const ParticleP2D &x1, const ParticleP2D &x2 )
{
    // Construct trajectory interpolation
    double dt = x2[0] - x1[0];
    TrajectoryRep1D trep[2];
    trep[0].construct( dt, x1[1], x1[2], x2[1], x2[2] );
    trep[1].construct( dt, x1[3], x1[4], x2[3], x2[4] );

    //std::cout << "trep[0]:\n";
    //trep[0].debug_print();
    //std::cout << "trep[1]:\n";
    //trep[1].debug_print();

    // Solve for intersections
    double K[3];
    //std::cout << "Solving for trep[" << crd << "] = " << val << ":\n";
    int nroots = trep[crd].solve( K, val );
    //std::cout << "found " << nroots << " roots: ";

    // Save intersection points
    double x[2], v[2];
    for( int b = 0; b < nroots; b++ ) {
	//std::cout << K[b] << " ";
	trep[0].coord( x[0], v[0], K[b] );
	trep[1].coord( x[1], v[1], K[b] );
	intsc.push_back( ParticleP2D( x1[0]+dt*K[b], x[0], v[0], x[1], v[1] ) );
    }
    //std::cout << "\n";

    return( nroots );
}


int ParticlePCyl::get_derivatives( double t, const double *x, double *dxdt, void *data )
{
    Vec3D E, B, xc( x[0], x[2], 0.0 );
    ParticleIteratorData *pi = (ParticleIteratorData *)data;

#ifdef DEBUG_PARTICLE_DERIVATIVES
    std::cout << "Particle get_derivatives query\n";
    std::cout << "  t    = " << t << "\n";    
    std::cout << "  x    = " << std::setw(12) << x[0] << " " 
	      << std::setw(12) << x[1] << " " 
	      << std::setw(12) << x[2] << " " 
	      << std::setw(12) << x[3] << " " 
	      << std::setw(12) << x[4] << "\n";
#endif

    if( x[2] <= 0.0 ) {
#ifdef DEBUG_PARTICLE_DERIVATIVES
	std::cout << "  Breaking iteration step\n";
#endif
	return( IBSIMU_DERIV_ERROR );
    }
    
    if( pi->_efield )
	E = (*pi->_efield)( xc );
    if( pi->_bfield )
	B = (*pi->_bfield)( xc );

#ifdef DEBUG_PARTICLE_DERIVATIVES
    std::cout << "  E    = " << E << "\n";
    std::cout << "  B    = " << B << "\n";
#endif
    
    /* Positions: dx/dt = vx, dr/dt = vr */
    dxdt[0] = x[1];
    dxdt[2] = x[3];
    
    /* Velocities:
     * dvx/dt = ax
     * dvr/dt = ar+r*(dtheta/dt)^2 
     * d^2theta/dt^2 = (a_theta-dr/dt*dtheta/dt)/r 
     */
    dxdt[1] = pi->_qm * (E[0] + x[3]*B[2] - x[2]*x[4]*B[1]);
    dxdt[3] = pi->_qm * (E[1] + x[2]*x[4]*B[0] - x[1]*B[2]) + x[2]*x[4]*x[4];
    if( x[2] == 0.0 )
        dxdt[4] = 0.0;
    else
        dxdt[4] = (pi->_qm * (x[1]*B[1] - x[3]*B[0]) - 2.0*x[3]*x[4]) / x[2];
    
#ifdef DEBUG_PARTICLE_DERIVATIVES
    std::cout << "  dxdt = " 
	      << std::setw(12) << dxdt[0] << " "
	      << std::setw(12) << dxdt[1] << " "
	      << std::setw(12) << dxdt[2] << " "
	      << std::setw(12) << dxdt[3] << " "
	      << std::setw(12) << dxdt[4] << "\n";
#endif

    return( GSL_SUCCESS );
}


int ParticlePCyl::trajectory_intersections_at_plane( std::vector<ParticlePCyl> &intsc, 
						     int crd, double val,
						     const ParticlePCyl &x1, const ParticlePCyl &x2 )
{
    // Construct trajectory interpolation
    double dt = x2[0] - x1[0];
    TrajectoryRep1D trep[2];
    trep[0].construct( dt, x1[1], x1[2], x2[1], x2[2] );
    trep[1].construct( dt, x1[3], x1[4], x2[3], x2[4] );

    // Solve for intersections
    double K[3];
    int nroots = trep[crd].solve( K, val );
    
    // Save intersection points
    double x[2], v[2];
    for( int b = 0; b < nroots; b++ ) {
	trep[0].coord( x[0], v[0], K[b] );
	trep[1].coord( x[1], v[1], K[b] );
	intsc.push_back( ParticlePCyl( x1[0]+K[b]*dt, x[0], v[0], x[1], v[1], x1[5]+K[b]*(x2[5]-x1[5]) ) );
    }

    return( nroots );
}


int ParticleP3D::get_derivatives( double t, const double *x, double *dxdt, void *data )
{
    Vec3D E, B, xc( x[0], x[2], x[4] );
    ParticleIteratorData *pi = (ParticleIteratorData *)data;

    if( pi->_efield )
	E = (*pi->_efield)( xc );
    if( pi->_bfield )
	B = (*pi->_bfield)( xc );

    /*
    std::cout << "\nget_derivatives query:\n";
    std::cout << "t=" << t << "\n";
    std::cout << "x=(" << xc << ")\n";
    std::cout << "E=(" << E << ")\n";
    std::cout << "B=(" << B << ")\n";
    std::cout << "qm= " << pi->_qm << "\n";    
    */

    /* Positions: dx/dt = vx, dy/dt = vy, dz/dt = vz */
    dxdt[0] = x[1];
    dxdt[2] = x[3];
    dxdt[4] = x[5];
    
    /* Velocities:
     * dvx/dt = ax
     * dvy/dt = ay
     * dvz/dt = az
     */
    dxdt[1] = pi->_qm * (E[0] + x[3]*B[2] - x[5]*B[1]);
    dxdt[3] = pi->_qm * (E[1] + x[5]*B[0] - x[1]*B[2]);
    dxdt[5] = pi->_qm * (E[2] + x[1]*B[1] - x[3]*B[0]);
    
    //std::cout << "dxdt=(" 
    //<< dxdt[0] << " "
    //<< dxdt[1] << " "
    //<< dxdt[2] << " "
    //<< dxdt[3] << " "
    //<< dxdt[4] << ")\n";

    return( GSL_SUCCESS );
}


int ParticleP3D::trajectory_intersections_at_plane( std::vector<ParticleP3D> &intsc, 
						    int crd, double val,
						    const ParticleP3D &x1, const ParticleP3D &x2 )
{
    // Construct trajectory interpolation
    double dt = x2[0] - x1[0];
    TrajectoryRep1D trep[3];
    trep[0].construct( dt, x1[1], x1[2], x2[1], x2[2] );
    trep[1].construct( dt, x1[3], x1[4], x2[3], x2[4] );
    trep[2].construct( dt, x1[5], x1[6], x2[5], x2[6] );

    // Solve for intersections
    double K[3];
    int nroots = trep[crd].solve( K, val );
    
    // Save intersection points
    double x[3], v[3];
    for( int b = 0; b < nroots; b++ ) {
	trep[0].coord( x[0], v[0], K[b] );
	trep[1].coord( x[1], v[1], K[b] );
	trep[2].coord( x[2], v[2], K[b] );
	intsc.push_back( ParticleP3D( x1[0]+dt*K[b], x[0], v[0], x[1], v[1], x[2], v[2] ) );
    }

    return( nroots );
}

















