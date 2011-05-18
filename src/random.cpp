/*! \file qrandom.cpp
 *  \brief Source code for qrandom.cpp
 */

/* Copyright (c) 2005-2009 Taneli Kalvas. All rights reserved.
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

#include <cmath>
#include "qrandom.hpp"
#include "error.hpp"


const double __rgauss_const_A = 1.0/sqrt(2.0);
const double __rgauss_const_B = 1.0/sqrt(2.0*M_PI);


/* Subroutine providing cumulative distribution function of gaussian at x. */
double Random::rgauss_f_func( double x, void *params )
{
    double *R = (double *)params;
    return( 0.5*( 1.0+erf(__rgauss_const_A) ) - *R );
}

/* Subroutine providing derivative of cumulative distribution function
 * of gaussian at x. */
double Random::rgauss_df_func( double x, void *params )
{
    return( __rgauss_const_B*exp( -x*x*__rgauss_const_A ) );
}

/* Subroutine providing cumulative distribution function of gaussian
 * and its derivative at x. */
void Random::rgauss_fdf_func( double x, void *params, double *f, double *df )
{
    double *R = (double *)params;

    *f  = 0.5*( 1.0+erf(x*__rgauss_const_A) ) - *R;
    *df = __rgauss_const_B*exp( -0.5*x*x );
}


/* Subroutine to filter uniformly distributed R to gaussian
 * distributed variable x by iterative solver. */
double Random::transform_gaussian( double R )
{
    int iter = 0;
    double x, x0;

    x = 0.0;
    _fdf.params = (void *)&R;
    gsl_root_fdfsolver_set( _solver, &_fdf, x );
    do {
        gsl_root_fdfsolver_iterate( _solver );
        x0 = x;
        x = gsl_root_fdfsolver_root( _solver );
        if( fabs(x - x0) < 1e-6 )
	    break;
	iter++;
    } while( iter < 100 );
    if( iter == 100 )
        throw( Error( ERROR_LOCATION, "too many iterations" ) );
    return( x );
}


Random::Random( int n )
    : _N(n)
{
    _fdf.f   = &rgauss_f_func;
    _fdf.df  = &rgauss_df_func;
    _fdf.fdf = &rgauss_fdf_func;
    _solver  = gsl_root_fdfsolver_alloc( gsl_root_fdfsolver_newton );
}


Random::~Random()
{
    gsl_root_fdfsolver_free( _solver );
}





/* ********************************************************************************* */



QRandom::QRandom( int n )
    : Random(n)
{
    _qrng = gsl_qrng_alloc( gsl_qrng_sobol, _N );
}


QRandom::~QRandom()
{
    gsl_qrng_free( _qrng );
}


void QRandom::get_gaussian( double *x )
{
    int i;

    gsl_qrng_get( _qrng, x );
    for( i = 0; i < _N; i++ )
	x[i] = transform_gaussian( x[i] );
}


void QRandom::get_part_gaussian( bool *gaussian, double *x )
{
    int i;

    gsl_qrng_get( _qrng, x );
    for( i = 0; i < _N; i++ ) {
	if( gaussian[i] )
	    x[i] = transform_gaussian( x[i] );
    }
}


void QRandom::get( double *x )
{
    gsl_qrng_get( _qrng, x );
}




/* ********************************************************************************* */



MTRandom::MTRandom( int n )
    : Random(n)
{
    _rng = gsl_rng_alloc( gsl_rng_mt19937 );
}


MTRandom::~MTRandom()
{
    gsl_rng_free( _rng );
}


void MTRandom::get_gaussian( double *x )
{
    for( int i = 0; i < _N; i++ )
	x[i] = gsl_rng_uniform( _rng );

    for( int i = 0; i < _N; i++ )
	x[i] = transform_gaussian( x[i] );
}


void MTRandom::get_part_gaussian( bool *gaussian, double *x )
{
    for( int i = 0; i < _N; i++ )
	x[i] = gsl_rng_uniform( _rng );

    for( int i = 0; i < _N; i++ ) {
	if( gaussian[i] )
	    x[i] = transform_gaussian( x[i] );
    }
}


void MTRandom::get( double *x )
{
    for( int i = 0; i < _N; i++ )
	x[i] = gsl_rng_uniform( _rng );
}



















