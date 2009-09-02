#include "gauss_seidel.hpp"
#include <cmath>
#include <limits>


bool gauss_seidel( const CRowMatrix &mat, const Vector &rhs, Vector &sol,
		   uint32_t &imax, double &eps, double w )
{
    uint32_t i;       /* Number of iterations */
    uint32_t j, k, l; /* Row, data and column pointers */
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
















