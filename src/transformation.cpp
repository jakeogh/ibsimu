/*! \file transformation.cpp
 *  \brief Affine transformation
 */

/* Copyright (c) 2010 Taneli Kalvas. All rights reserved.
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


#include "transformation.hpp"


double Transformation::determinant( void ) const
{
    return( + x[3] * x[6] * x[9] * x[12]
            - x[2] * x[7] * x[9] * x[12]
            - x[3] * x[5] * x[10] * x[12]
            + x[1] * x[7] * x[10] * x[12]
            + x[2] * x[5] * x[11] * x[12]
            - x[1] * x[6] * x[11] * x[12]
            - x[3] * x[6] * x[8] * x[13]
            + x[2] * x[7] * x[8] * x[13]
            + x[3] * x[4] * x[10] * x[13]
            - x[0] * x[7] * x[10] * x[13]
            - x[2] * x[4] * x[11] * x[13]
            + x[0] * x[6] * x[11] * x[13]
            + x[3] * x[5] * x[8] * x[14]
            - x[1] * x[7] * x[8] * x[14]
            - x[3] * x[4] * x[9] * x[14]
            + x[0] * x[7] * x[9] * x[14]
            + x[1] * x[4] * x[11] * x[14]
            - x[0] * x[5] * x[11] * x[14]
            - x[2] * x[5] * x[8] * x[15]
            + x[1] * x[6] * x[8] * x[15]
            + x[2] * x[4] * x[9] * x[15]
            - x[0] * x[6] * x[9] * x[15]
            - x[1] * x[4] * x[10] * x[15]
            + x[0] * x[5] * x[10] * x[15] );
}


Transformation Transformation::inverse( void ) const
{
    double idet = 1.0/determinant();
    Transformation ret( x[6]*x[11]*x[13] - x[7]*x[10]*x[13] + x[7]*x[9]*x[14] - 
			x[5]*x[11]*x[14] - x[6]*x[9]*x[15] + x[5]*x[10]*x[15],
			x[3]*x[10]*x[13] - x[2]*x[11]*x[13] - x[3]*x[9]*x[14] + 
			x[1]*x[11]*x[14] + x[2]*x[9]*x[15] - x[1]*x[10]*x[15],
			x[2]*x[7]*x[13] - x[3]*x[6]*x[13] + x[3]*x[5]*x[14] - 
			x[1]*x[7]*x[14] - x[2]*x[5]*x[15] + x[1]*x[6]*x[15],
			x[3]*x[6]*x[9] - x[2]*x[7]*x[9] - x[3]*x[5]*x[10] + 
			x[1]*x[7]*x[10] + x[2]*x[5]*x[11] - x[1]*x[6]*x[11],
			x[7]*x[10]*x[12] - x[6]*x[11]*x[12] - x[7]*x[8]*x[14] + 
			x[4]*x[11]*x[14] + x[6]*x[8]*x[15] - x[4]*x[10]*x[15],
			x[2]*x[11]*x[12] - x[3]*x[10]*x[12] + x[3]*x[8]*x[14] - 
			x[0]*x[11]*x[14] - x[2]*x[8]*x[15] + x[0]*x[10]*x[15],
			x[3]*x[6]*x[12] - x[2]*x[7]*x[12] - x[3]*x[4]*x[14] + 
			x[0]*x[7]*x[14] + x[2]*x[4]*x[15] - x[0]*x[6]*x[15],
			x[2]*x[7]*x[8] - x[3]*x[6]*x[8] + x[3]*x[4]*x[10] - 
			x[0]*x[7]*x[10] - x[2]*x[4]*x[11] + x[0]*x[6]*x[11],
			x[5]*x[11]*x[12] - x[7]*x[9]*x[12] + x[7]*x[8]*x[13] - 
			x[4]*x[11]*x[13] - x[5]*x[8]*x[15] + x[4]*x[9]*x[15],
			x[3]*x[9]*x[12] - x[1]*x[11]*x[12] - x[3]*x[8]*x[13] + 
			x[0]*x[11]*x[13] + x[1]*x[8]*x[15] - x[0]*x[9]*x[15],
			x[1]*x[7]*x[12] - x[3]*x[5]*x[12] + x[3]*x[4]*x[13] - 
			x[0]*x[7]*x[13] - x[1]*x[4]*x[15] + x[0]*x[5]*x[15],
			x[3]*x[5]*x[8] - x[1]*x[7]*x[8] - x[3]*x[4]*x[9] + 
			x[0]*x[7]*x[9] + x[1]*x[4]*x[11] - x[0]*x[5]*x[11],
			x[6]*x[9]*x[12] - x[5]*x[10]*x[12] - x[6]*x[8]*x[13] + 
			x[4]*x[10]*x[13] + x[5]*x[8]*x[14] - x[4]*x[9]*x[14],
			x[1]*x[10]*x[12] - x[2]*x[9]*x[12] + x[2]*x[8]*x[13] - 
			x[0]*x[10]*x[13] - x[1]*x[8]*x[14] + x[0]*x[9]*x[14],
			x[2]*x[5]*x[12] - x[1]*x[6]*x[12] - x[2]*x[4]*x[13] + 
			x[0]*x[6]*x[13] + x[1]*x[4]*x[14] - x[0]*x[5]*x[14],
			x[1]*x[6]*x[8] - x[2]*x[5]*x[8] + x[2]*x[4]*x[9] - 
			x[0]*x[6]*x[9] - x[1]*x[4]*x[10] + x[0]*x[5]*x[10] );
    ret *= idet;
    return( ret );
}


const Transformation &Transformation::operator*=( double s )
{
    for( size_t i = 0; i < 16; i++ )
        x[i] *= s;
    return( *this );
}


Transformation Transformation::operator*( const Transformation &m ) const
{
    Transformation r;

    r[ 0] = x[ 0]*m[ 0] + x[ 1]*m[ 4] + x[ 2]*m[ 8] + x[ 3]*m[12];
    r[ 1] = x[ 0]*m[ 1] + x[ 1]*m[ 5] + x[ 2]*m[ 9] + x[ 3]*m[13];
    r[ 2] = x[ 0]*m[ 2] + x[ 1]*m[ 6] + x[ 2]*m[10] + x[ 3]*m[14];
    r[ 3] = x[ 0]*m[ 3] + x[ 1]*m[ 7] + x[ 2]*m[11] + x[ 3]*m[15];

    r[ 4] = x[ 4]*m[ 0] + x[ 5]*m[ 4] + x[ 6]*m[ 8] + x[ 7]*m[12];
    r[ 5] = x[ 4]*m[ 1] + x[ 5]*m[ 5] + x[ 6]*m[ 9] + x[ 7]*m[13];
    r[ 6] = x[ 4]*m[ 2] + x[ 5]*m[ 6] + x[ 6]*m[10] + x[ 7]*m[14];
    r[ 7] = x[ 4]*m[ 3] + x[ 5]*m[ 7] + x[ 6]*m[11] + x[ 7]*m[15];

    r[ 8] = x[ 8]*m[ 0] + x[ 9]*m[ 4] + x[10]*m[ 8] + x[11]*m[12];
    r[ 9] = x[ 8]*m[ 1] + x[ 9]*m[ 5] + x[10]*m[ 9] + x[11]*m[13];
    r[10] = x[ 8]*m[ 2] + x[ 9]*m[ 6] + x[10]*m[10] + x[11]*m[14];
    r[11] = x[ 8]*m[ 3] + x[ 9]*m[ 7] + x[10]*m[11] + x[11]*m[15];

    r[12] = x[12]*m[ 0] + x[13]*m[ 4] + x[14]*m[ 8] + x[15]*m[12];
    r[13] = x[12]*m[ 1] + x[13]*m[ 5] + x[14]*m[ 9] + x[15]*m[13];
    r[14] = x[12]*m[ 2] + x[13]*m[ 6] + x[14]*m[10] + x[15]*m[14];
    r[15] = x[12]*m[ 3] + x[13]*m[ 7] + x[14]*m[11] + x[15]*m[15];

    return( r );
}


Vec4D Transformation::operator*( const Vec4D &v ) const
{
    Vec4D r;

    r[0] =  x[0]*v[0] +  x[1]*v[1] +  x[2]*v[2] +  x[3]*v[3];
    r[1] =  x[4]*v[0] +  x[5]*v[1] +  x[6]*v[2] +  x[7]*v[3];
    r[2] =  x[8]*v[0] +  x[9]*v[1] + x[10]*v[2] + x[11]*v[3];
    r[3] = x[12]*v[0] + x[13]*v[1] + x[14]*v[2] + x[15]*v[3];

    return( r );
}


Vec4D Transformation::operator%( const Vec4D &v ) const
{
    Vec4D r;

    r[0] =  x[0]*v[0] +  x[4]*v[1] +  x[8]*v[2] + x[12]*v[3];
    r[1] =  x[1]*v[0] +  x[5]*v[1] +  x[9]*v[2] + x[13]*v[3];
    r[2] =  x[2]*v[0] +  x[6]*v[1] + x[10]*v[2] + x[14]*v[3];
    r[3] =  x[3]*v[0] +  x[7]*v[1] + x[11]*v[2] + x[15]*v[3];

    return( r );
}


Vec3D Transformation::transform_point( const Vec3D &xin ) const
{
    Vec4D r = *this * Vec4D( xin[0], xin[1], xin[2], 1.0 );

    return( Vec3D( r[0], r[1], r[2] ) );
}


Vec3D Transformation::inv_transform_point( const Vec3D &xin ) const
{
    Transformation t = this->inverse();
    Vec4D r = t * Vec4D( xin[0], xin[1], xin[2], 1.0 );

    return( Vec3D( r[0], r[1], r[2] ) );
}


Vec3D Transformation::transform_vector( const Vec3D &xin ) const
{
    Vec4D r = *this * Vec4D( xin[0], xin[1], xin[2], 0.0 );

    return( Vec3D( r[0], r[1], r[2] ) );
}


Vec3D Transformation::inv_transform_vector( const Vec3D &xin ) const
{
    Transformation t = this->inverse();
    Vec4D r = t * Vec4D( xin[0], xin[1], xin[2], 0.0 );

    return( Vec3D( r[0], r[1], r[2] ) );
}
