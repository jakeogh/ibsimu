/*! \file axisymmetricvectorfield.cpp
 *  \brief Axisymmetric magnetic field
 */

/* Copyright (c) 2012 Taneli Kalvas. All rights reserved.
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


#include "axisymmetricvectorfield.hpp"


AxisymmetricVectorField::AxisymmetricVectorField( geom_mode_e geom_mode, 
						  std::vector<double> z,
						  std::vector<double> Bz )
{
    if( geom_mode != MODE_3D && geom_mode != MODE_CYL )
	throw( Error( ERROR_LOCATION, "unsupported geometry mode" ) );
    _geom_mode = geom_mode;
    if( z.size() != Bz.size() )
	throw( Error( ERROR_LOCATION, "non-matching Bz and z arrays" ) );
    _spline = gsl_spline_alloc( gsl_interp_cspline, z.size() );
    _accel = gsl_interp_accel_alloc();
    gsl_spline_init( _spline, &z[0], &Bz[0], z.size() );
}


AxisymmetricVectorField::AxisymmetricVectorField( const AxisymmetricVectorField &f )
{
    _geom_mode = f._geom_mode;
    _spline = gsl_spline_alloc( gsl_interp_cspline, f._spline->size );
    _accel = gsl_interp_accel_alloc();
    gsl_spline_init( _spline, f._spline->x, f._spline->y, f._spline->size );
}


AxisymmetricVectorField &AxisymmetricVectorField::operator=( const AxisymmetricVectorField &f )
{
    gsl_spline_free( _spline );
    gsl_interp_accel_free( _accel );

    _geom_mode = f._geom_mode;
    _spline = gsl_spline_alloc( gsl_interp_cspline, f._spline->size );
    _accel = gsl_interp_accel_alloc();
    gsl_spline_init( _spline, f._spline->x, f._spline->y, f._spline->size );
    return( *this );
}


AxisymmetricVectorField::~AxisymmetricVectorField()
{
    gsl_spline_free( _spline );
    gsl_interp_accel_free( _accel );
}


const Vec3D AxisymmetricVectorField::operator()( const Vec3D &x ) const
{
    Vec3D B;

    if( _geom_mode == MODE_3D ) {
	double Bz0  = gsl_spline_eval( _spline, x[2], _accel );
	double Bzp  = gsl_spline_eval_deriv( _spline, x[2], _accel );
	double Bzpp = gsl_spline_eval_deriv2( _spline, x[2], _accel );

	double rr = x[0]*x[0] + x[1]*x[1];
	double Bz = Bz0 - rr*Bzpp/4.0;
	double Br = -sqrt(rr)*Bzp/2.0;
	double theta = atan2( x[1], x[0] );
	B = Vec3D( Br*cos(theta), Br*sin(theta), Bz );
    } else /* _geom_mode == MODE_CYL */ {
	double Bz0  = gsl_spline_eval( _spline, x[0], _accel );
	double Bzp  = gsl_spline_eval_deriv( _spline, x[0], _accel );
	double Bzpp = gsl_spline_eval_deriv2( _spline, x[0], _accel );

	double Bz = Bz0 - x[1]*x[1]*Bzpp/4.0;
	double Br = -x[1]*Bzp/2.0;
	B = Vec3D( Bz, Br, 0.0 );
    }

    return( B );
}


