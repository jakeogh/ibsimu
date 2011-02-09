/*! \file bicgstab.cpp
 *  \brief Stabilized Biconjugate Gradient solver
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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

#include <iostream>
#include <iomanip>
#include <sstream>

#include "compmath.hpp"
#include "bicgstab.hpp"
#include "ibsimu.hpp"
#include "statusprint.hpp"


void bicgstab( const Matrix &mat, const Vector &rhs, Vector &sol,
	       const Precond &pc, uint32_t &imax, double &eps )
{
    // Checks
    if( mat.columns() != mat.rows() )
	throw( ErrorDim( ERROR_LOCATION, "matrix not square" ) );
    if( mat.rows() != rhs.size() )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    double resid, omega = 0, alpha = 0, beta, rho_1, rho_2 = 0;
    Vector p, phat, s, shat, t, v;
    double norm_rhs = norm2(rhs);
    if( sol.size() != mat.columns() ) {
	sol.resize( mat.columns() );
	sol.clear();
    }

    Vector r = mat * sol;
    r = rhs - r;
    Vector rtilde = r;

    if( norm_rhs == 0.0 )
	norm_rhs = 1;
  
    if( (resid = norm2(r) / norm_rhs) <= eps ) {
	eps = resid;
	imax = 0;
	return;
    }

    StatusPrint sp;
    if( ibsimu.get_verbose_output() ) {
	std::stringstream ss;
	ss << "  " << std::setw(5) << 0 << " " << std::scientific << std::setw(20) << resid;
	sp.print( ss.str() );
    }

    uint32_t i;
    for( i = 1; i <= imax; i++ ) {
	rho_1 = dot_prod( rtilde, r );
	if( rho_1 == 0 ) {
	    eps = norm2(r) / norm_rhs;
	    imax = i;
	    break;
	}
	if( i == 1 )
	    p = r;
	else {
	    beta = (rho_1/rho_2) * (alpha/omega);
	    p = r + beta * (p - omega * v);
	}
	pc.solve( phat, p );
	v = mat * phat;
	alpha = rho_1 / dot_prod( rtilde, v );
	s = r - alpha * v;
	if( (resid = norm2(s)/norm_rhs) < eps ) {
	    sol += alpha * phat;
	    eps = resid;
	    imax = i;
	    break;
	}
	pc.solve( shat, s );
	t = mat * shat;
	omega = dot_prod( t, s ) / dot_prod( t, t );
	sol += alpha * phat + omega * shat;
	r = s - omega * t;

	rho_2 = rho_1;
	if( (resid = norm2(r) / norm_rhs) < eps ) {
	    eps = resid;
	    imax = i;
	    break;
	}
	if( comp_isnan( resid ) || omega == 0 ) {
	    throw( Error( ERROR_LOCATION, "convergence failure" ) );
	}

	if( ibsimu.get_verbose_output() ) {
	    std::stringstream ss;
	    ss << "  " << std::setw(5) << i << " " << std::scientific << std::setw(20) << resid;
	    sp.print( ss.str() );
	}
    }

    if( i > imax ) {
	eps = resid;
	imax = i;
    }
    
    return;
}

