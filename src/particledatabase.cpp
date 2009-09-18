/*! \file particledatabase.cpp
 *  \brief Source code for particledatabase.cpp
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

#include "particledatabase.hpp"
#include "qrandom.hpp"


void ParticleDataBase2D::add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
						    double v, double dvp, double dvt, 
						    double x1, double y1, double x2, double y2 )
{
    if( verbose_output )
	std::cout << "Defining a 2d beam\n";

    _particles.reserve( _particles.size()+N );

    QRandom qrng( 2 );

    m *= MASS_U;
    q *= CHARGE_E;
    double s = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    _rhosum += J/v;
    double IQ = J*s/N; // IQ is in units A/m for 2D
    double vt[2], dv[2];
    ParticleP2D x;
    x[0] = 0.0;
    Vec3D transverse( x2-x1, y2-y1, 0.0 );
    transverse /= transverse.norm2();
    Vec3D parallel( transverse[1], -transverse[0], 0.0 );

    double Isum = 0.0;
    for( uint32_t a = 0; a < N; a++ ) {
	x[1] = x1 + (x2-x1)*(a+0.5)/((double)N);
	x[3] = y1 + (y2-y1)*(a+0.5)/((double)N);

	qrng.get_gaussian( vt );
	dv[0] = transverse[0]*dvt*vt[0] + parallel[0]*dvp*vt[1];
	dv[1] = transverse[1]*dvt*vt[0] + parallel[1]*dvp*vt[1];
	x[2] = parallel[0]*v + dv[0];
	x[4] = parallel[1]*v + dv[1];

	_particles.push_back( Particle2D( IQ, q, m, x ) );
	Isum += IQ;
    }

    if( verbose_output )
	std::cout << "  Total beam current " << Isum << " A/m\n";
}


void ParticleDataBase2D::add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
						  double E, double Tp, double Tt, 
						  double x1, double y1, double x2, double y2 )
{
    add_2d_beam_with_velocity( N, J, q, m, sqrt(2.0*E*CHARGE_E/(m*MASS_U)), 
			       sqrt(Tp*CHARGE_E/(m*MASS_U)), 
			       sqrt(Tt*CHARGE_E/(m*MASS_U)), 
			       x1, y1, x2, y2 );
}


void ParticleDataBaseCyl::add_2d_beam_with_velocity( uint32_t N, double J, double q, double m, 
						     double v, double dvp, double dvt, 
						     double x1, double y1, double x2, double y2 )
{
    if( verbose_output )
	std::cout << "Defining a cylindrical beam\n";

    _particles.reserve( _particles.size()+N );

    QRandom qrng( 3 );

    m *= MASS_U;
    q *= CHARGE_E;
    double s = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    _rhosum += J/v;
    double IQ = 2.0*M_PI*J*s/N;
    double vt[3], dv[3];
    ParticlePCyl x;
    x[0] = 0.0;
    Vec3D transverse( x2-x1, y2-y1, 0.0 );
    transverse /= transverse.norm2();
    Vec3D parallel( transverse[1], -transverse[0], 0.0 );

    double Isum = 0.0;
    for( uint32_t a = 0; a < N; a++ ) {
	x[1] = x1 + (x2-x1)*(a+0.5)/((double)N);
	x[3] = y1 + (y2-y1)*(a+0.5)/((double)N);

	qrng.get_gaussian( vt );
	dv[0] = transverse[0]*dvt*vt[0] + parallel[0]*dvp*vt[1];
	dv[1] = transverse[1]*dvt*vt[0] + parallel[1]*dvp*vt[1];
	dv[2] = dvt*vt[2];
	x[2] = parallel[0]*fabs(v) + dv[0];
	x[4] = parallel[1]*fabs(v) + dv[1];
	if( x[3] == 0.0 )
	    x[5] = 0.0;
	else
	    x[5] = dv[2]/x[3];

	_particles.push_back( ParticleCyl( IQ*x[3], q, m, x ) );
	Isum += IQ*x[3];
    }

    if( verbose_output )
	std::cout << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBaseCyl::add_2d_beam_with_energy( uint32_t N, double J, double q, double m, 
						   double E, double Tp, double Tt, 
						   double x1, double y1, double x2, double y2 )
{
    add_2d_beam_with_velocity( N, J, q, m, sqrt(2.0*E*CHARGE_E/(m*MASS_U)), 
			       sqrt(Tp*CHARGE_E/(m*MASS_U)), 
			       sqrt(Tt*CHARGE_E/(m*MASS_U)), 
			       x1, y1, x2, y2 );
}




void ParticleDataBase3D::add_cylindrical_beam_with_velocity( uint32_t N, double J, double q, double m, 
							     double v, double dvp, double dvt, Vec3D c, 
							     Vec3D dir1, Vec3D dir2, double r )
{
    if( verbose_output )
	std::cout << "Defining a cylindrical beam\n";

    _particles.reserve( _particles.size()+N );

    // Random number generator for two positions and three velocities (gaussian)
    QRandom qrng( 5 );
    bool randmask[5] = {false, false, true, true, true};
    double qx[5];
    double px[6];

    m *= MASS_U;
    q *= CHARGE_E;
    _rhosum += J/v;
    double IQ = M_PI*r*r*J/N;

    // Calculate and check base vectors
    Vec3D dir3 = cross( dir1, dir2 );
    dir2 = cross( dir1, dir3 );
    dir1.normalize();
    dir2.normalize();
    dir3.normalize();
    if( dir1[0] != dir1[0] || dir2[0] != dir2[0] || dir3[0] != dir3[0] ) {
	throw( Error( ERROR_LOCATION, "invalid direction vectors" ) );
    }

    // Prepare particle
    ParticleP3D x;
    x[0] = 0.0;

    double Isum = 0.0;
    uint32_t a = 0;
    while( a < N ) {

	qrng.get_part_gaussian( randmask, qx );

	// Calculate in natural (dir1,dir2,dir3) coordinates
	px[0] = -r + 2.0*r*qx[0];
	px[1] = -r + 2.0*r*qx[1];
	px[2] = 0.0;
	px[3] = dvt*qx[2];
	px[4] = dvt*qx[3];
	px[5] = dvp*qx[4] + v;

	// Reject points outside unit radius circle
	if( px[0]*px[0] + px[1]*px[1] > r*r )
	    continue;

	// Map to world coordinates
	x[1] = dir1[0]*px[0] + dir2[0]*px[1] + dir3[0]*px[2] + c[0];
	x[2] = dir1[0]*px[3] + dir2[0]*px[4] + dir3[0]*px[5];
	x[3] = dir1[1]*px[0] + dir2[1]*px[1] + dir3[1]*px[2] + c[1];
	x[4] = dir1[1]*px[3] + dir2[1]*px[4] + dir3[1]*px[5];
	x[5] = dir1[2]*px[0] + dir2[2]*px[1] + dir3[2]*px[2] + c[2];
	x[6] = dir1[2]*px[3] + dir2[2]*px[4] + dir3[2]*px[5];

	_particles.push_back( Particle3D( IQ, q, m, x ) );
	Isum += IQ;
	a++;
    }

    if( verbose_output )
	std::cout << "  Total beam current " << Isum << " A\n";
}


void ParticleDataBase3D::add_cylindrical_beam_with_energy( uint32_t N, double J, double q, double m, 
							   double E, double Tp, double Tt, Vec3D c,
							   Vec3D dir1, Vec3D dir2, double r )
{
    add_cylindrical_beam_with_velocity( N, J, q, m, 
					sqrt(2.0*E*CHARGE_E/(m*MASS_U)), 
					sqrt(Tp*CHARGE_E/(m*MASS_U)),
					sqrt(Tt*CHARGE_E/(m*MASS_U)),
					c, dir1, dir2, r );
}















