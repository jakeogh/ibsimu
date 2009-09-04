#include "bicgstab_solver.hpp"
#include "bicgstab.hpp"
#include "ilu0_precond.hpp"
#include "verbose.hpp"
#include "timer.hpp"
#include "error.hpp"


BiCGSTABSolver::BiCGSTABSolver( double eps, int imax,
				double newton_Reps, double newton_dXeps, int newton_imax )
    : _eps(eps), _imax(imax), _newton_Reps(newton_Reps), _newton_dXeps(newton_dXeps), _newton_imax(newton_imax)
{
    if( _imax <= 0 || _newton_imax <= 0 )
	throw( ErrorDim( ERROR_LOCATION, "invalid number of iterations" ) );
    if( eps <= 0.0 || newton_Reps <= 0.0 || newton_dXeps <= 0.0 )
	throw( ErrorDim( ERROR_LOCATION, "invalid accuracy request" ) );
}


void BiCGSTABSolver::solve( const Problem &p, Vector &X ) const
{
    int imax;
    double eps;
    Timer t;
    if( p.linear() ) {

	// Linear solver
	if( verbose_output )
	    std::cout << "  Using ILU0-BiCGSTAB solver\n";

	const Matrix *A;
	const Vector *B;
	p.get_vecmat( &A, &B );
	ILU0_Precond pc( *A );
	imax = _imax;
	eps = _eps;
	bicgstab( *A, *B, X, pc, imax, eps );

	if( verbose_output ) {
	    std::cout << "  iterations = " << imax << " (max " << _imax << ")\n";
	    std::cout << "  eps = " << eps << " (requested " << _eps << ")\n";
	}

    } else {

	// Nonlinear solver (Newton-Raphson)
	if( verbose_output )
	    std::cout << "  Using Newton-Raphson ILU0-BiCGSTAB solver\n";

	int32_t a;
	int imax_sum = 0;
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

	    if( verbose_output )
		std::cout << "    " 
			  << std::setw(5) << a << " " 
			  << std::setw(14) << accX << " " 
			  << std::setw(14) << accR << "\n";

	    if( accR < _newton_Reps || accX < _newton_dXeps || imax_sum >= _imax )
		break;
	}

	if( verbose_output ) {
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
    if( verbose_output )
	std::cout << "  time used = " << t << "\n";
}















