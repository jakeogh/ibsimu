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


void UMFPACKSolver::umfpack_solve( CColMatrix &mat, const Vector &rhs, Vector &sol )
{
    int     status;
    void   *Symbolic;
    void   *Numeric;

    // Sort
    //mat.order_ascending();
    if( !mat.check_ascending() )
	throw( Error( ERROR_LOCATION, "matrix not in ascending order" ) );

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

	CColMatrix colmatrix( *A );
	umfpack_solve( colmatrix, *B, X );

        if( verbose_output )
            std::cout << "  Done\n";

    } else {

        // Nonlinear solver (Newton-Raphson)
        if( verbose_output )
            std::cout << "  Using Newton-Raphson UMFPACK solver\n";

	std::cout << "Not implemented yet\n";
	exit( 1 );
	
    }

    t.stop();
    if( verbose_output )
        std::cout << "  time used = " << t << "\n";
}
