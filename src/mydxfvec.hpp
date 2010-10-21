/*! \file mydxfvec.hpp
 *  \brief Header file for mydxfvec.hpp
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

#ifndef MYDXF_VEC_HPP
#define MYDXF_VEC_HPP 1


#include <math.h>
#include <iostream>
#include <iomanip>
#include "error.hpp"


/*! \brief DXF vector class.
 *
 *  A basic three dimensional vector for DXF operations.
 */
class MyDXFVec
{
    double p[3];

public:

    MyDXFVec( double x = 0.0, double y = 0.0, double z = 0.0 ) {
	p[0] = x;
	p[1] = y;
	p[2] = z;
    }
    ~MyDXFVec() {}

    double &operator[]( int i ) { return( p[i] ); }
    const double &operator[]( int i ) const { return( p[i] ); }
    double &operator()( int i ) { return( p[i] ); }
    const double &operator()( int i ) const { return( p[i] ); }

    /*! \brief %Vector addition
     */
    MyDXFVec operator+( const MyDXFVec &vec ) const { 
        return( MyDXFVec( p[0] + vec[0], 
                       p[1] + vec[1],
                       p[2] + vec[2] ) );
    }

    /*! \brief %Vector difference
     */
    MyDXFVec operator-( const MyDXFVec &vec ) const {
        return( MyDXFVec( p[0] - vec[0],
                       p[1] - vec[1],
                       p[2] - vec[2] ) );
    } 

    /*! \brief %Vector accumulation
     */
    MyDXFVec &operator+=( const MyDXFVec &vec ) { 
        p[0] += vec[0];
        p[1] += vec[1];
        p[2] += vec[2];
        return( *this );
    }

    /*! \brief Dot product
     */
    double operator*( const MyDXFVec &vec ) const { 
        return( p[0] * vec[0] +
                p[1] * vec[1] +
                p[2] * vec[2] );
    }

    MyDXFVec operator*( double x ) { 
        return( MyDXFVec( x*p[0], x*p[1], x*p[2] ) );
    }

    MyDXFVec &operator*=( double x ) { 
        p[0] *= x;
        p[1] *= x;
        p[2] *= x;
        return( *this );
    }

    MyDXFVec &operator/=( double x ) { 
        double div = 1.0/x;
        p[0] *= div;
        p[1] *= div;
        p[2] *= div;
        return( *this );
    }

    bool operator!=( const MyDXFVec &x ) { 
        if( p[0] != x.p[0] || p[1] != x.p[1] || p[2] != x.p[2] )
            return( true );
        return( false ); 
    }

    bool operator==( const MyDXFVec &x ) { 
        if( p[0] == x.p[0] && p[1] == x.p[1] && p[2] == x.p[2] )
            return( true );
        return( false ); 
    }

    MyDXFVec &operator=( const MyDXFVec &x ) { 
        p[0] = x[0];
        p[1] = x[1];
        p[2] = x[2];
        return( *this );
    }

    MyDXFVec &operator=( const double &x ) { 
        p[0] = x;
        p[1] = x;
        p[2] = x;
        return( *this );
    }

    /*! \brief Normalize vector
     */
    void normalize() {
        double inv_norm = 1.0/sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
        p[0] *= inv_norm;
        p[1] *= inv_norm;
        p[2] *= inv_norm;
    }

    /*! \brief Returns 2-norm of vector
     *
     *  \f$ ||x||_2 = \sqrt{ \Sigma_{i=1}^n x_i^2 } \f$
     */
    double norm2() const {
        return( sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] ) );
    }

    /*! \brief Returns square of 2-norm of vector
     *
     *  \f$ (||x||_2)^2 = \Sigma_{i=1}^n x_i^2 \f$
     */
    double ssqr() const {
        return( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
    }

    /*! \brief Cross product
     */
    friend MyDXFVec cross( const MyDXFVec &vec1, const MyDXFVec &vec2 );

    friend double norm2( const MyDXFVec &vec );

    friend MyDXFVec operator*( double x, const MyDXFVec &vec );
    friend std::ostream &operator<<( std::ostream &os, const MyDXFVec &vec );
};

inline double norm2( const MyDXFVec &vec ) {
    return( vec.norm2() );
}

inline MyDXFVec cross( const MyDXFVec &vec1, const MyDXFVec &vec2 ) { 
    return( MyDXFVec( vec1[1] * vec2[2] - vec1[2] * vec2[1], 
                   vec1[2] * vec2[0] - vec1[0] * vec2[2],
                   vec1[0] * vec2[1] - vec1[1] * vec2[0] ) );
}


inline MyDXFVec operator*( double x, const MyDXFVec &vec )
{
    return( MyDXFVec( x*vec.p[0], x*vec.p[1], x*vec.p[2] ) );
}


inline std::ostream &operator<<( std::ostream &os, const MyDXFVec &vec ) 
{
    os << std::setw(12) << to_string(vec[0]).substr(0,12) << " ";
    os << std::setw(12) << to_string(vec[1]).substr(0,12) << " ";
    os << std::setw(12) << to_string(vec[2]).substr(0,12);
    return( os );
}


#endif
