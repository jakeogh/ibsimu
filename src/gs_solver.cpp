/*! \file gs_solver.cpp
 *  \brief Source code for gs_solver.cpp
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

#include <limits>
#include "gs_solver.hpp"
#include "ibsimu.hpp"
#include "timer.hpp"
#include "error.hpp"


GSSolver::GSSolver( double w,
		    double eps, uint32_t imax,
		    double newton_Reps, 
		    double newton_dXeps, 
		    uint32_t newton_imax )
    : _eps(eps), _imax(imax), _w(w), _newton_Reps(newton_Reps), 
      _newton_dXeps(newton_dXeps), _newton_imax(newton_imax)
{
    
}



void GSSolver::solve( const Problem &p, Vector &X )
{
    uint32_t imax;
    double eps;
    Timer t;
    if( p.linear() ) {

        // Linear solver
        if( ibsimu.get_verbose_output() )
            std::cout << "  Using Gauss-Seidel solver\n";

        const Matrix *A; // Matrix needs to be modified for smooth edges
        const Vector *B;
        p.get_vecmat( &A, &B );
        imax = _imax;
        eps = _eps;
        gauss_seidel( *A, *B, X, imax, eps, _w );

        if( ibsimu.get_verbose_output() ) {
            std::cout << "  iterations = " << imax << " (max " << _imax << ")\n";
            std::cout << "  eps = " << eps << " (requested " << _eps << ")\n";
        }
    } else {
	throw( ErrorUnimplemented( ERROR_LOCATION, "non-linear solver unimplemented" ) );
    }

    t.stop();
    if( ibsimu.get_verbose_output() )
	std::cout << "  time used = " << t << "\n";
}


void GSSolver::reset( void )
{

}


bool GSSolver::gauss_seidel( const CRowMatrix &mat, const Vector &rhs, Vector &sol,
			     uint32_t &imax, double &eps, double w )
{
    uint32_t i;          /* Number of iterations */
    int j, k, l;         /* Row, data and column pointers */
    double dia = 0;      /* Diagonal element */
    double x;            /* Temporary variable */
    double d;            /* Change */
    double d_max;        /* Maximum change */
    double w2 = 1.0-w;

    sol.resize( rhs.size() );

    i = 0;
    d_max = std::numeric_limits<double>::infinity();
    while( i < imax && d_max > eps ) {
	d_max = 0.0;
	/* Go through the rows */
	for( j = 0; j < mat.rows(); j++ ) {
	    /* For each row the new Xi value is Xi = (B - sum(Aij*Xj) ) / Aii,
	     * where i != j */
	    x = 0.0;
	    /* Go through non-zero cells */
	    for( k = mat.ptr(j); k < mat.ptr(j+1); k++ ) {
		/* Save diagonal element */
		l = mat.col(k);
		if( l == j )
		    dia = mat.val(k);
		else
 		    x += mat.val(k) * sol(l);
	    }
	    x = w2*sol(j) - w*(x-rhs(j)) / dia;
	    if( !finite( x ) ) {
		imax = i;
		eps = x;
		return( false );
	    }
	    d = fabs( x - sol(j) );
	    sol(j) = x;
	    if( d > d_max )
		d_max = d;
	}

	i++;
    }

    imax = i;
    eps = d_max;
    return( true );
}







