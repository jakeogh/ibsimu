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
















