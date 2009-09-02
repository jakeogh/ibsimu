#include <cmath>
#include "qrandom.hpp"
#include "error.hpp"


const double _QRandom_A = 1.0/sqrt(2.0);
const double _QRandom_B = 1.0/sqrt(2.0*M_PI);


/* Subroutine providing cumulative distribution function of gaussian at x. */
double _QRandom_f_func( double x, void *params )
{
    double *R = (double *)params;
    return( 0.5*( 1.0+erf(_QRandom_A) ) - *R );
}

/* Subroutine providing derivative of cumulative distribution function
 * of gaussian at x. */
double _QRandom_df_func( double x, void *params )
{
    return( _QRandom_B*exp( -x*x*_QRandom_A ) );
}

/* Subroutine providing cumulative distribution function of gaussian
 * and its derivative at x. */
void _QRandom_fdf_func( double x, void *params, double *f, double *df )
{
    double *R = (double *)params;

    *f  = 0.5*( 1.0+erf(x*_QRandom_A) ) - *R;
    *df = _QRandom_B*exp( -0.5*x*x );
}


/* Subroutine to filter uniformly distributed R to gaussian
 * distributed variable x by iterative solver. */
double QRandom::transform_gaussian( double R )
{
    int iter = 0;
    double x, x0;

    x = 0.0;
    fdf.params = (void *)&R;
    gsl_root_fdfsolver_set( solver, &fdf, x );
    do {
        gsl_root_fdfsolver_iterate( solver );
        x0 = x;
        x = gsl_root_fdfsolver_root( solver );
        if( fabs(x - x0) < 1e-6 )
	    break;
	iter++;
    } while( iter < 100 );
    if( iter == 100 )
        throw( Error( ERROR_LOCATION, "too many iterations" ) );
    return( x );
}


QRandom::QRandom( int n )
{
    N = n;
    qrng = gsl_qrng_alloc( gsl_qrng_sobol, N );
    fdf.f = &_QRandom_f_func;
    fdf.df = &_QRandom_df_func;
    fdf.fdf = &_QRandom_fdf_func;
    solver = gsl_root_fdfsolver_alloc( gsl_root_fdfsolver_newton );
}


QRandom::~QRandom()
{
    gsl_root_fdfsolver_free( solver );
    gsl_qrng_free( qrng );
}


void QRandom::get_gaussian( double *x )
{
    int i;

    gsl_qrng_get( qrng, x );
    for( i = 0; i < N; i++ )
	x[i] = transform_gaussian( x[i] );
}


void QRandom::get_part_gaussian( bool *gaussian, double *x )
{
    int i;

    gsl_qrng_get( qrng, x );
    for( i = 0; i < N; i++ ) {
	if( gaussian[i] )
	    x[i] = transform_gaussian( x[i] );
    }
}


void QRandom::get( double *x )
{
    gsl_qrng_get( qrng, x );
}














