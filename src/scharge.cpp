/*! \file scharge.cpp
 *  \brief Space charge deposition functions
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

#include "scharge.hpp"
#include "ibsimu.hpp"


//#define DEBUG_SCHARGE 1


void scharge_finalize( ScalarField &scharge )
{
    if( ibsimu.get_verbose_output() )
	ibsimu.vout() << "  Finalizing space charge density map\n";

    switch( scharge.geom_mode() ) {
    case MODE_2D:
    {
	// Convert charge map to space charge density map
	scharge /= (scharge.h()*scharge.h());

	// Correct boundaries
	for( int32_t i = 0; i < scharge.size(0); i++ ) {
	    scharge( i, 0 ) *= 2.0;
	    scharge( i, scharge.size(1)-1 ) *= 2.0;
	}
	for( int32_t j = 0; j < scharge.size(1); j++ ) {
	    scharge( 0, j ) *= 2.0;
	    scharge( scharge.size(0)-1, j ) *= 2.0;
	}
	break;
    }
    case MODE_CYL:
    {
	// Convert charge map to space charge density map
 	for( int32_t i = 0; i < scharge.size(0); i++ ) {
	    for( int32_t j = 0; j < scharge.size(1); j++ ) {
		if( j == 0 ) {
                    double rj2 = scharge.h()+scharge.origo(1);
                    scharge( i, j ) /= (M_PI*scharge.h()*(rj2*rj2));
                } else {
                    double rj1 = (j-0.5)*scharge.h()+scharge.origo(1);
                    double rj2 = (j+0.5)*scharge.h()+scharge.origo(1);
                    scharge( i, j ) /= (M_PI*scharge.h()*(rj2*rj2-rj1*rj1));
                }
	    }
	}
	// Correct boundaries
	for( int32_t i = 0; i < scharge.size(0); i++ ) {
	    scharge( i, 0 ) *= 2.0;
	    scharge( i, scharge.size(1)-1 ) *= 2.0;
	}
	for( int32_t j = 0; j < scharge.size(1); j++ ) {
	    scharge( 0, j ) *= 2.0;
	    scharge( scharge.size(0)-1, j ) *= 2.0;
	}
	break;
    }
    case MODE_3D:
    {
	// Convert charge map to space charge density map
	scharge /= (scharge.h()*scharge.h()*scharge.h());

	// Correct boundaries
 	for( int32_t i = 0; i < scharge.size(0); i++ ) {
	    for( int32_t j = 0; j < scharge.size(1); j++ ) {
		scharge( i, j, 0 ) *= 2.0;
		scharge( i, j, scharge.size(2)-1 ) *= 2.0;
	    }
	}
	for( int32_t i = 0; i < scharge.size(0); i++ ) {
	    for( int32_t k = 0; k < scharge.size(2); k++ ) {
		scharge( i, 0, k ) *= 2.0;
		scharge( i, scharge.size(1)-1, k ) *= 2.0;
	    }
	}	
	for( int32_t j = 0; j < scharge.size(1); j++ ) {
	    for( int32_t k = 0; k < scharge.size(2); k++ ) {
		scharge( 0, j, k ) *= 2.0;
		scharge( scharge.size(0)-1, j, k ) *= 2.0;
	    }
	}	
	break;
    }
    default:
    {
	throw( Error( ERROR_LOCATION, "unsupported dimension number" ) );
    }
    }
}


void scharge_add_from_trajectory( ScalarField &scharge, pthread_mutex_t *mutex, double IQ, 
				  const ParticleP2D &x1, const ParticleP2D &x2 )
{
    double x[2];
    double t[2];
    int i[2];

#ifdef DEBUG_SCHARGE
    std::cout << "Calculating space charge\n";
    std::cout << "x1 = " << x1 << "\n";
    std::cout << "x2 = " << x2 << "\n";
#endif
    for( size_t a = 0; a < 2; a++ ) {
	x[a] = 0.5*( x1[2*a+1] + x2[2*a+1] );
	i[a] = (int)floor( ( x[a]-scharge.origo(a) ) * scharge.div_h() );
	t[a] = ( x[a]-(i[a]*scharge.h()+scharge.origo(a)) ) * scharge.div_h();
	
#ifdef DEBUG_SCHARGE
	std::cout << "a = " << a << "\n";
	std::cout << "x = " << x[a] << "\n";
	std::cout << "i = " << i[a] << "\n";
	std::cout << "t = " << t[a] << "\n";
#endif

	// Add charge to boundaries when over simulation area
	if( i[a] < 0 ) {
	    i[a] = 0;
	    t[a] = 0.0;
	} else if( i[a] >= scharge.size(a)-1 ) {
	    i[a] = scharge.size(a)-2;
	    t[a] = 1.0;
	}
    }

    /*
    if( i[0] < 0 || i[0] >= scharge.size(0)-1 ||
	i[1] < 0 || i[1] >= scharge.size(1)-1 ) {
	throw( Error( ERROR_LOCATION, "coordinates out of range at ("
		      + to_string(i[0]) + ", "
		      + to_string(i[1]) + ")" ) );
    }
    */

    double Q = IQ*(x2[0]-x1[0]); // Q = I*dt
#ifdef DEBUG_SCHARGE
    std::cout << "IQ = " << IQ << "\n";
    std::cout << "dt = " << (x2[0]-x1[0]) << "\n";
    std::cout << "Q = " << Q << "\n\n";
#endif
    int p = scharge.size(0)*i[1] + i[0];
    pthread_mutex_lock( mutex );
    scharge( p )                   += (1.0-t[0])*(1.0-t[1])*Q;
    scharge( p+scharge.size(0) )   += (1.0-t[0])*t[1]*Q;
    scharge( p+1 )                 += t[0]*(1.0-t[1])*Q;
    scharge( p+1+scharge.size(0) ) += t[0]*t[1]*Q;
    pthread_mutex_unlock( mutex );
}


void scharge_add_from_trajectory( ScalarField &scharge, pthread_mutex_t *mutex, double IQ, 
				  const ParticlePCyl &x1, const ParticlePCyl &x2 )
{
    double x[2];
    double t[2];
    int i[2];

    // x-direction
    x[0] = 0.5*( x1[1] + x2[1] );
    i[0] = (int)floor( ( x[0]-scharge.origo(0) ) * scharge.div_h() );
    t[0] = ( x[0]-(i[0]*scharge.h()+scharge.origo(0)) ) * scharge.div_h();

    // r-direction
    x[1] = 0.5*( x1[3] + x2[3] );
    i[1] = (int)floor( ( x[1]-scharge.origo(1) ) * scharge.div_h() );
    double rj1 = i[1]*scharge.h()+scharge.origo(1);
    double rj2 = rj1+scharge.h();
    rj1 = rj1*rj1;
    rj2 = rj2*rj2;
    t[1] = (x[1]*x[1]-rj1) / (rj2-rj1);

    for( size_t a = 0; a < 2; a++ ) {
	// Add charge to boundaries when over simulation area
	if( i[a] < 0 ) {
	    i[a] = 0;
	    t[a] = 0.0;
	} else if( i[a] >= scharge.size(a)-1 ) {
	    i[a] = scharge.size(a)-2;
	    t[a] = 1.0;
	}
    }

    double Q = IQ*(x2[0]-x1[0]); // Q = I*dt
    int p = scharge.size(0)*i[1] + i[0];
    pthread_mutex_lock( mutex );
    scharge( p )                   += (1.0-t[0])*(1.0-t[1])*Q;
    scharge( p+scharge.size(0) )   += (1.0-t[0])*t[1]*Q;
    scharge( p+1 )                 += t[0]*(1.0-t[1])*Q;
    scharge( p+1+scharge.size(0) ) += t[0]*t[1]*Q;    
    pthread_mutex_unlock( mutex );
}


void scharge_add_from_trajectory( ScalarField &scharge, pthread_mutex_t *mutex, double IQ, 
				  const ParticleP3D &x1, const ParticleP3D &x2 )
{
    double x[3];
    double t[3];
    int i[3];

    for( size_t a = 0; a < 3; a++ ) {
	x[a] = 0.5*( x1[2*a+1] + x2[2*a+1] );
	i[a] = (int)floor( ( x[a]-scharge.origo(a) ) * scharge.div_h() );
	t[a] = ( x[a]-(i[a]*scharge.h()+scharge.origo(a)) ) * scharge.div_h();

	// Add charge to boundaries when over simulation area
	if( i[a] < 0 ) {
	    i[a] = 0;
	    t[a] = 0.0;
	} else if( i[a] >= scharge.size(a)-1 ) {
	    i[a] = scharge.size(a)-2;
	    t[a] = 1.0;
	}
    }

    double Q = IQ*(x2[0]-x1[0]); // Q = I*dt
    int p = scharge.size(0)*scharge.size(1)*i[2] + scharge.size(0)*i[1] + i[0];

    pthread_mutex_lock( mutex );
    scharge( p )                   += (1.0-t[0])*(1.0-t[1])*(1.0-t[2])*Q;
    scharge( p+scharge.size(0) )   += (1.0-t[0])*t[1]*(1.0-t[2])*Q;
    scharge( p+1 )                 += t[0]*(1.0-t[1])*(1.0-t[2])*Q;
    scharge( p+1+scharge.size(0) ) += t[0]*t[1]*(1.0-t[2])*Q;

    p += scharge.size(0)*scharge.size(1);
    scharge( p )                   += (1.0-t[0])*(1.0-t[1])*t[2]*Q;
    scharge( p+scharge.size(0) )   += (1.0-t[0])*t[1]*t[2]*Q;
    scharge( p+1 )                 += t[0]*(1.0-t[1])*t[2]*Q;
    scharge( p+1+scharge.size(0) ) += t[0]*t[1]*t[2]*Q;
    pthread_mutex_unlock( mutex );
}



















