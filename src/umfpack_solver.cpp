#include <cstdlib>
#include <iostream>
#include <umfpack.h>
#include "umfpack_solver.hpp"
#include "timer.hpp"
#include "verbose.hpp"
#include "error.hpp"



UMFPACKSolver::UMFPACKSolver( double newton_Reps, double newton_dXeps, uint32_t newton_imax )
    : _newton_Reps(newton_Reps), _newton_dXeps(newton_dXeps), _newton_imax(newton_imax)
{

}


void UMFPACKSolver::umfpack_error( const std::string func, int status )
{
    if( status == UMFPACK_ERROR_n_nonpositive )
	throw( Error( ERROR_LOCATION, func + ": n less than or equal to zero" ) );
    else if( status == UMFPACK_ERROR_invalid_matrix )
	throw( Error( ERROR_LOCATION, func + ": invalid matrix" ) );
    else if( status == UMFPACK_ERROR_out_of_memory )
	throw( Error( ERROR_LOCATION, func + ": memory allocation error" ) );
    else if( status == UMFPACK_ERROR_argument_missing )
	throw( Error( ERROR_LOCATION, func + ": argument missing" ) );
    else if( status == UMFPACK_ERROR_internal_error )
	throw( Error( ERROR_LOCATION, func + ": internal error" ) );
    else if( status == UMFPACK_WARNING_singular_matrix )
	throw( Error( ERROR_LOCATION, func + ": singular matrix" ) );
    else if( status == UMFPACK_ERROR_invalid_Symbolic_object )
	throw( Error( ERROR_LOCATION, func + ": invalid symbolic object" ) );
    else if( status == UMFPACK_ERROR_invalid_system )
	throw( Error( ERROR_LOCATION, func + ": invalid system" ) );
    else if( status == UMFPACK_ERROR_different_pattern )
	throw( Error( ERROR_LOCATION, func + ": different pattern" ) );
    else if( status == UMFPACK_ERROR_invalid_Numeric_object )
	throw( Error( ERROR_LOCATION, func + ": invalid numeric object" ) );
    else if( status != UMFPACK_OK )
	throw( Error( ERROR_LOCATION, "unknown error in " + func ) );
}


void UMFPACKSolver::umfpack_solve( const CColMatrix &mat, const Vector &rhs, Vector &sol )
{
    int     status;
    void   *Symbolic;
    void   *Numeric;

    sol.resize( rhs.size() );

    // Sort
    //mat.order_ascending();
    //if( !mat.check_ascending() )
    //	throw( Error( ERROR_LOCATION, "matrix not in ascending order" ) );

    // Make decomposition
    status = umfpack_di_symbolic( mat.columns(), mat.rows(), 
				  &mat.ptr(0), &mat.row(0), &mat.val(0), 
				  &Symbolic, (double *)NULL, (double *)NULL );
    if( status != UMFPACK_OK )
	umfpack_error( "umfpack_di_symbolic", status );

    status = umfpack_di_numeric( &mat.ptr(0), &mat.row(0), &mat.val(0), 
				 Symbolic, &Numeric, (double *)NULL, 
				 (double *)NULL );
    if( status != UMFPACK_OK )
	umfpack_error( "umfpack_di_numeric", status );

    // Solve using decomposition
    status = umfpack_di_solve( UMFPACK_A, &mat.ptr(0), &mat.row(0), &mat.val(0), 
			       sol.get_data(), rhs.get_data(), Numeric, 
			       (double *)NULL, (double *)NULL );
    if( status != UMFPACK_OK )
	umfpack_error( "umfpack_di_solve", status );
    
    umfpack_di_free_symbolic( &Symbolic );
    umfpack_di_free_numeric( &Numeric );
}


void UMFPACKSolver::solve( const Problem &p, Vector &X ) const
{
    Timer t;
    if( p.linear() ) {

        // Linear solver
        if( verbose_output )
            std::cout << "  Using UMFPACK solver\n";

	const Matrix *A;
        const Vector *B;
        p.get_vecmat( &A, &B );

	CColMatrix Acol( *A );
	umfpack_solve( Acol, *B, X );

        if( verbose_output )
            std::cout << "  Done\n";

    } else {

        // Nonlinear solver (Newton-Raphson)
        if( verbose_output )
            std::cout << "  Using Newton-Raphson UMFPACK solver\n";

	int32_t a;
	const Matrix *J;
	const Vector *R;
	double accR = 0.0, accX = 0.0;
	Vector dX;

	if( verbose_output )
	    std::cout << "    " 
		      << std::setw(5) << "Iter" << " " 
		      << std::setw(14) << "Step size" << " " 
		      << std::setw(14) << "Residual" << "\n";

	for( a = 0; a < _newton_imax; a++ ) {
	    // Calculate dX = J^{-1}*R	    
	    p.get_resjac( &J, &R, X );
	    CColMatrix Jcol( *J );
	    dX.clear();
	    umfpack_solve( Jcol, *R, dX );

	    // Take step
	    X -= dX;

	    // Check for convergence
	    accR = max_abs( *R );
	    accX = max_abs( dX );

	    if( verbose_output )
		std::cout << "    " 
			  << std::setw(5) << a << " " 
			  << std::setw(14) << accX << " " 
			  << std::setw(14) << accR << "\n";

	    if( accR < _newton_Reps || accX < _newton_dXeps )
		break;
	}

	if( verbose_output ) {
	    if( accR < _newton_Reps || accX < _newton_dXeps )
		std::cout << "  Newton-Raphson converged\n";
	    else
		std::cout << "  Maximum number of Newton-Raphson iterations\n";
	}
    }

    t.stop();
    if( verbose_output )
        std::cout << "  time used = " << t << "\n";
}
