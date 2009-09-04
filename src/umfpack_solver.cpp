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


void UMFPACKSolver::umfpack_solve( CColMatrix &mat, const Vector &rhs, Vector &sol )
{
    int     status;
    void   *Symbolic;
    void   *Numeric;

    // Sort
    mat.order_ascending();

    // Make decomposition
    status = umfpack_di_symbolic( mat.columns(), mat.rows(), 
				  &mat.ptr(0), &mat.row(0), &mat.val(0), 
				  &Symbolic, (double *)NULL, (double *)NULL );
    if( status != UMFPACK_OK )
	throw( Error( ERROR_LOCATION, "error in umfpack_di_symbolic" ) );
    status = umfpack_di_numeric( &mat.ptr(0), &mat.row(0), &mat.val(0), 
				 Symbolic, &Numeric, (double *)NULL, 
				 (double *)NULL );
    if( status == UMFPACK_WARNING_singular_matrix )
	throw( Error( ERROR_LOCATION, "singular matrix" ) );
    else if( status != UMFPACK_OK )
	throw( Error( ERROR_LOCATION, "error in umfpack_di_numeric" ) );

    // Solve using decomposition
    status = umfpack_di_solve( UMFPACK_A, &mat.ptr(0), &mat.row(0), &mat.val(0), 
			       sol.get_data(), rhs.get_data(), Numeric, 
			       (double *)NULL, (double *)NULL );
    if( status == UMFPACK_ERROR_out_of_memory )
	throw( Error( ERROR_LOCATION, "memory allocation error" ) );
    else if( status == UMFPACK_WARNING_singular_matrix )
	throw( Error( ERROR_LOCATION, "singular matrix" ) );
    else if( status != UMFPACK_OK )
	throw( Error( ERROR_LOCATION, "error in umfpack_di_solve" ) );
    
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
