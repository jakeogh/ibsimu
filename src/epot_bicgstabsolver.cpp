/*! \file epot_bicgstabsolver.cpp
 *  \brief BiCGSTAB matrix solver for electric potential problem
 */

/* Copyright (c) 2005-2012 Taneli Kalvas. All rights reserved.
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
#include "bicgstab.hpp"
#include "hbio.hpp"
#include "ilu0_precond.hpp"
#include "epot_bicgstabsolver.hpp"
#include "ibsimu.hpp"


EpotBiCGSTABSolver::EpotBiCGSTABSolver( Geometry &geom, 
					double eps, 
					uint32_t imax,
					double newton_Reps, 
					double newton_dXeps, 
					uint32_t newton_imax,
					bool gnewton )
    : EpotMatrixSolver(geom), _eps(eps), _imax(imax), _iter(0), _res(0.0), _gnewton(gnewton),
      _newton_Reps(newton_Reps), _newton_dXeps(newton_dXeps), _newton_imax(newton_imax)
{
    if( eps <= 0.0 || newton_Reps <= 0.0 || newton_dXeps <= 0.0 )
        throw( ErrorDim( ERROR_LOCATION, "invalid accuracy request" ) );

    _pc = new ILU0_Precond;
}


EpotBiCGSTABSolver::EpotBiCGSTABSolver( Geometry &geom, std::istream &s )
    : EpotMatrixSolver(geom,s)
{
    _pc = new ILU0_Precond;
    throw( ErrorUnimplemented( ERROR_LOCATION ) );

}


void EpotBiCGSTABSolver::save( std::ostream &s ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


EpotBiCGSTABSolver::~EpotBiCGSTABSolver()
{
    if( _pc )
	delete _pc;
}


void EpotBiCGSTABSolver::set_preconditioner( Precond &pc )
{
    if( _pc )
	delete _pc;
    _pc = pc.copy();
}


void EpotBiCGSTABSolver::set_gnewton( bool enable ) 
{
    _gnewton = enable;
}


void EpotBiCGSTABSolver::set_eps( double eps ) 
{
    if( eps <= 0.0 )
        throw( ErrorDim( ERROR_LOCATION, "invalid accuracy request" ) );
    _eps = eps;
}


void EpotBiCGSTABSolver::set_imax( uint32_t imax ) 
{
    _imax = imax;
}


void EpotBiCGSTABSolver::set_newton_imax( uint32_t newton_imax ) 
{
    _newton_imax = newton_imax;
}


void EpotBiCGSTABSolver::set_newton_residual_eps( double newton_Reps ) 
{
    if( newton_Reps <= 0.0 )
        throw( ErrorDim( ERROR_LOCATION, "invalid accuracy request" ) );
    _newton_Reps = newton_Reps;
}


void EpotBiCGSTABSolver::set_newton_step_eps( double newton_dXeps ) 
{
    if( newton_dXeps <= 0.0 )
        throw( ErrorDim( ERROR_LOCATION, "invalid accuracy request" ) );
    _newton_dXeps = newton_dXeps;
}


double EpotBiCGSTABSolver::get_residual( void ) const
{
    return( _res );
}


uint32_t EpotBiCGSTABSolver::get_iter( void ) const
{
    return( _iter );
}


void EpotBiCGSTABSolver::reset_problem( void )
{
    reset_matrix();
    _pc->clear();
}


void EpotBiCGSTABSolver::subsolve( MeshScalarField &epot, const MeshScalarField &scharge )
{
    uint32_t imax;
    double eps;

    // Preprocess and set starting guess
    preprocess( epot, scharge );
    Vector X;
    set_initial_guess( epot, X );

    if( linear() ) {

	ibsimu.message( 1 ) << "Using BiCGSTAB-" << _pc->typestring() << " solver("
			    << " imax = " << _imax
			    << ", eps = " << _eps
			    << " )\n";
	ibsimu.flush();
	    
	// Fetch matrix form of problem
	const CRowMatrix *A;
	const Vector *B;
	get_vecmat( &A, &B );

	if( !_pc->is_prepared() )
	    _pc->prepare( *A );
	_pc->construct( *A );

        imax = _imax;
        eps = _eps / _res_coef;
        bicgstab( *A, *B, X, *_pc, imax, eps );
	_iter = imax;
	_res = _res_coef * eps;

	if( _iter == _imax )
	    ibsimu.message( 1 ) << "Maximum number of iteration rounds done.\n";
	ibsimu.message( 1 ) << "residual error = " << _res << "\n";
	ibsimu.message( 1 ) << "iterations = " << _iter << "\n";
	ibsimu.flush();

    } else {

	int32_t a;
	uint32_t imax_sum = 0;
        const CRowMatrix *J;
        const Vector *R;
        double accR = 0.0, accX = 0.0;
        Vector dX;

	if( _gnewton ) {

	    ibsimu.message( 1 ) << "Using Globally convergent Newton-Raphson BiCGSTAB-" 
				<< _pc->typestring() << " solver("
				<< " imax = " << _imax
				<< ", eps = " << _eps
				<< ", newton_imax = " << _newton_imax
				<< ", newton_reps = " << _newton_Reps
				<< ", newton_dxeps = " << _newton_dXeps
				<< " )\n";
	    ibsimu.message( 1 ) << "  " 
				<< std::setw(5)  << "Round" << " " 
				<< std::setw(8)  << "Iter" << " " 
				<< std::setw(14) << "Step size" << " " 
				<< std::setw(14) << "Step fac" << " " 
				<< std::setw(14) << "Residual" << "\n";
	    ibsimu.flush();
	    
	    // Globally convergent Newton-Raphson
            Vector Xold( X.size() );

            // First jacobian and residual
            get_resjac( &J, &R, X );
            double f = ssqr( *R );

	    if( !_pc->is_prepared() )
		_pc->prepare( *J );

	    for( a = 0; a < (int)_newton_imax; a++ ) {

                // Calculate dX = J^{-1}*R
		_pc->construct( *J );
                imax = _imax - imax_sum;
                eps = _eps / _res_coef;
                dX.clear();
                bicgstab( *J, *R, dX, *_pc, imax, eps );
                imax_sum += imax;

                // Search for acceptable step for which residual decreases
                double t = 2.0;
                double fold = f;
                Xold = X;
                while( f >= fold ) {

                    t *= 0.5;
                    X = Xold - t*dX;
                    get_resjac( &J, &R, X );
                    f = ssqr( *R );
                    if( t <= std::numeric_limits<double>::epsilon() )
                        break;
                }

                // Check for convergence
                accR = max_abs( *R );
                accX = t*max_abs( dX );

		ibsimu.message( 1 ) << "  " 
				    << std::setw(5)  << a << " " 
				    << std::setw(8)  << imax << " " 
				    << std::setw(14) << accX << " " 
				    << std::setw(14) << t << " " 
				    << std::setw(14) << accR << "\n";
		ibsimu.flush();
                
                if( accR < _newton_Reps || (t == 1.0 && accX < _newton_dXeps) || imax_sum >= _imax )
                    break;
            }

	} else {

	    ibsimu.message( 1 ) << "Using Newton-Raphson BiCGSTAB-" << _pc->typestring() << " solver("
				<< " imax = " << _imax
				<< ", eps = " << _eps
				<< ", newton_imax = " << _newton_imax
				<< ", newton_reps = " << _newton_Reps
				<< ", newton_dxeps = " << _newton_dXeps
				<< " )\n";
	    ibsimu.message( 1 ) << "  " 
				<< std::setw(5)  << "Round" << " " 
				<< std::setw(8)  << "Iter" << " " 
				<< std::setw(14) << "Step size" << " " 
				<< std::setw(14) << "Residual" << "\n";
	    ibsimu.flush();

	    for( a = 0; a < (int)_newton_imax; a++ ) {

		// Calculate dX = J^{-1}*R
		get_resjac( &J, &R, X );
		if( !_pc->is_prepared() )
		    _pc->prepare( *J );
		_pc->construct( *J );

		imax = _imax - imax_sum;
		eps = _eps / _res_coef;
		dX.clear();
		bicgstab( *J, *R, dX, *_pc, imax, eps );
		imax_sum += imax;
		
		// Take step
		X -= dX;

		// Check for convergence
		accR = max_abs( *R );
		accX = max_abs( dX );

		ibsimu.message( 1 ) << "  " 
				    << std::setw(5) << a << " " 
				    << std::setw(8) << imax << " " 
				    << std::setw(14) << accX << " " 
				    << std::setw(14) << _res_coef * accR << "\n";
		ibsimu.flush();
		
		if( accR < _newton_Reps || accX < _newton_dXeps || imax_sum >= _imax )
		    break;
	    }
        }

	_iter = imax_sum;
	_res = _res_coef * accR;

	if( accR < _newton_Reps || accX < _newton_dXeps )
	    ibsimu.message( 1 ) << "Newton-Raphson converged\n";
	else if( imax_sum >= _imax )
	    ibsimu.message( 1 ) << "Maximum number of BiCGSTAB iterations\n";
	else
	    ibsimu.message( 1 ) << "Maximum number of Newton-Raphson iterations\n";
	
	ibsimu.message( 1 ) << "residual error = " << _res << "\n";
	ibsimu.message( 1 ) << "total iterations = " << _iter << "\n";

    }

    // Postprocess and set solution
    set_solution( epot, X );
    postprocess();
}

    
void EpotBiCGSTABSolver::debug_print( std::ostream &os ) const
{
    EpotMatrixSolver::debug_print( os );
    os << "**EpotBiCGSTABSolver\n";


}
