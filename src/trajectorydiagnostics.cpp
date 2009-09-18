/*! \file trajectorydiagnostics.cpp
 *  \brief Source code for trajectorydiagnostics.cpp
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

#include <iostream>
#include <limits>
#include <cstring>
#include <cmath>
#include "trajectorydiagnostics.hpp"


void TrajectoryDiagnosticColumn::mirror( coordinate_axis_e axis, double level )
{
    size_t size = _data.size();
    _data.reserve( 2*size );

    // Handle X-axis
    if( axis == AXIS_X ) {
	if( _diag == DIAG_X ) {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( 2*level-_data[a] );
	} else if( _diag == DIAG_VX || _diag == DIAG_XP ) {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( -_data[a] );
	} else {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( _data[a] );
	}
    }

    // Handle Y-axis
    else if( axis == AXIS_Y || axis == AXIS_R ) {
	if( _diag == DIAG_Y || _diag == DIAG_R ) {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( 2*level-_data[a] );
	} else if( _diag == DIAG_VY || _diag == DIAG_VR || _diag == DIAG_YP ||_diag == DIAG_RP ) {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( -_data[a] );
	} else {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( _data[a] );
	}
    }

    // Handle Z-axis
    else if( axis == AXIS_Z ) {
	if( _diag == DIAG_Z ) {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( 2*level-_data[a] );
	} else if( _diag == DIAG_VZ || _diag == DIAG_ZP ) {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( -_data[a] );
	} else {
	    for( size_t a = 0; a < size; a++ )
		_data.push_back( _data[a] );
	}
    }
}


Emittance::Emittance( const std::vector<double> &x,
		      const std::vector<double> &xp,
		      const std::vector<double> &w )
{
    size_t N = x.size() < xp.size() ? 
	(x.size() < w.size() ? x.size() : w.size()) : 
	(xp.size() < w.size() ? xp.size() : w.size());

    // Calculate averages
    _wsum  = 0.0;
    _xave  = 0.0;
    _xpave = 0.0;
    for( size_t a = 0; a < N; a++ ) {
	_wsum  += w[a];
	_xave  += x[a]*w[a];
	_xpave += xp[a]*w[a];
    }
    _xave  = _xave  / _wsum;
    _xpave = _xpave / _wsum;

    // Calculate expectation values
    _x2  = 0.0;
    _xp2 = 0.0;
    _xxp = 0.0;
    for( size_t a = 0; a < N; a++ ) {
	_x2 += (x[a]-_xave)*(x[a]-_xave)*w[a];
	_xp2 += (xp[a]-_xpave)*(xp[a]-_xpave)*w[a];
	_xxp += (x[a]-_xave)*(xp[a]-_xpave)*w[a];
    }
    _x2  = _x2  / _wsum;
    _xp2 = _xp2 / _wsum;
    _xxp = _xxp / _wsum;

    // Calculate Twiss parameters
    _epsilon = sqrt( _xp2*_x2 - _xxp*_xxp );
    _alpha   = -_xxp/_epsilon;
    _beta    = _x2/_epsilon;
    _gamma   = _xp2/_epsilon;

    _angle = 0.5*atan( (2.0*_alpha) / (_beta - _gamma) );
    double H = 0.5*(_beta+_gamma);
    _rmajor = sqrt( 0.5*_epsilon ) * ( sqrt(H+1.0)+sqrt(H-1.0) );
    _rminor = sqrt( 0.5*_epsilon ) * ( sqrt(H+1.0)-sqrt(H-1.0) );
}














