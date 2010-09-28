/*! \file bicgstab_solver.cpp
 *  \brief Source code for bicgstab_solver.cpp
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

#include "bicgstab_solver.hpp"
#include "bicgstab.hpp"
#include "ilu0_precond.hpp"
#include "ibsimu.hpp"
#include "timer.hpp"
#include "error.hpp"


BiCGSTABSolver::BiCGSTABSolver( double eps, uint32_t imax,
				double newton_Reps, 
				double newton_dXeps, 
				uint32_t newton_imax )
    : _eps(eps), _imax(imax), _newton_Reps(newton_Reps), 
      _newton_dXeps(newton_dXeps), _newton_imax(newton_imax)
{
    if( _imax <= 0 || _newton_imax <= 0 )
	throw( ErrorDim( ERROR_LOCATION, "invalid number of iterations" ) );
    if( eps <= 0.0 || newton_Reps <= 0.0 || newton_dXeps <= 0.0 )
	throw( ErrorDim( ERROR_LOCATION, "invalid accuracy request" ) );
}


void BiCGSTABSolver::solve( const Problem &p, Vector &X )
{
    int imax;
    double eps;
    Timer t;
    if( p.linear() ) {

	// Linear solver
	if( ibsimu.get_verbose_output() )
	    std::cout << "  Using ILU0-BiCGSTAB solver\n";

	const Matrix *A;
	const Vector *B;
	p.get_vecmat( &A, &B );
	ILU0_Precond pc( *A );
	imax = _imax;
	eps = _eps;
	bicgstab( *A, *B, X, pc, imax, eps );

	if( ibsimu.get_verbose_output() ) {
	    std::cout << "  iterations = " << imax << " (max " << _imax << ")\n";
	    std::cout << "  eps = " << eps << " (requested " << _eps << ")\n";
	}

    } else {

	// Nonlinear solver (Newton-Raphson)
	if( ibsimu.get_verbose_output() )
	    std::cout << "  Using Newton-Raphson ILU0-BiCGSTAB solver\n";

	int32_t a;
	int imax_sum = 0;
	const Matrix *J;
	const Vector *R;
	double accR = 0.0, accX = 0.0;
	Vector dX;

	if( ibsimu.get_verbose_output() )
	    std::cout << "    " 
		      << std::setw(5) << "Iter" << " " 
		      << std::setw(14) << "Step size" << " " 
		      << std::setw(14) << "Residual" << "\n";

	for( a = 0; a < _newton_imax; a++ ) {
	    // Calculate dX = J^{-1}*R	    
	    p.get_resjac( &J, &R, X );
	    ILU0_Precond pc( *J );
	    imax = _imax - imax_sum;
	    eps = _eps;
	    dX.clear();
	    bicgstab( *J, *R, dX, pc, imax, eps );
	    imax_sum += imax;

	    // Take step
	    X -= dX;

	    // Check for convergence
	    accR = max_abs( *R );
	    accX = max_abs( dX );

	    if( ibsimu.get_verbose_output() )
		std::cout << "    " 
			  << std::setw(5) << a << " " 
			  << std::setw(14) << accX << " " 
			  << std::setw(14) << accR << "\n";

	    if( accR < _newton_Reps || accX < _newton_dXeps || imax_sum >= _imax )
		break;
	}

	if( ibsimu.get_verbose_output() ) {
	    if( accR < _newton_Reps || accX < _newton_dXeps )
		std::cout << "  Newton-Raphson converged\n";
	    else if( imax_sum >= _imax )
		std::cout << "  Maximum number of BiCGSTAB iterations\n";
	    else
		std::cout << "  Maximum number of Newton-Raphson iterations\n";

	    std::cout << "  total iterations = " << imax_sum << " (max " << _imax << ")\n";
	    std::cout << "  eps = " << eps << " (requested " << _eps << ")\n";
	}
    }

    t.stop();
    if( ibsimu.get_verbose_output() )
	std::cout << "  time used = " << t << "\n";
}


void BiCGSTABSolver::reset( void )
{

}














