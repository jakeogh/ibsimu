#include <limits>
#include "gs_solver.hpp"
#include "verbose.hpp"
#include "timer.hpp"
#include "error.hpp"


GSSolver::GSSolver( double w,
		    double eps, int imax,
		    double newton_Reps, 
		    double newton_dXeps, 
		    int newton_imax )
    : _eps(eps), _imax(imax), _w(w), _newton_Reps(newton_Reps), 
      _newton_dXeps(newton_dXeps), _newton_imax(newton_imax)
{
    
}



void GSSolver::solve( const Problem &p, Vector &X )
{
    int imax;
    double eps;
    Timer t;
    if( p.linear() ) {

        // Linear solver
        if( verbose_output )
            std::cout << "  Using Gauss-Seidel solver\n";

        const Matrix *A; // Matrix needs to be modified for smooth edges
        const Vector *B;
        p.get_vecmat( &A, &B );
        imax = _imax;
        eps = _eps;
        gauss_seidel( *A, *B, X, imax, eps, _w );

        if( verbose_output ) {
            std::cout << "  iterations = " << imax << " (max " << _imax << ")\n";
            std::cout << "  eps = " << eps << " (requested " << _eps << ")\n";
        }
    } else {
	throw( ErrorUnimplemented( ERROR_LOCATION, "non-linear solver unimplemented" ) );
    }

    t.stop();
    if( verbose_output )
	std::cout << "  time used = " << t << "\n";
}


void GSSolver::reset( void )
{

}


bool GSSolver::gauss_seidel( const CRowMatrix &mat, const Vector &rhs, Vector &sol,
			     int &imax, double &eps, double w )
{
    int i;          /* Number of iterations */
    int j, k, l;    /* Row, data and column pointers */
    double dia = 0; /* Diagonal element */
    double x;       /* Temporary variable */
    double d;       /* Change */
    double d_max;   /* Maximum change */
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


