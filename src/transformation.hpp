/*! \file transformation.hpp
 *  \brief Header file for affine transformation
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

#ifndef TRANSFORMATION_HPP
#define TRANSFORMATION_HPP 1


#include <string.h>
#include "vec3d.hpp"


/*! \brief Affine transformation for three dimensional space.
 *
 *  Affine transformation operating on Vec3D vectors and points.
 */
class Transformation
{

    double _x[16]; /* Matrix data in row first order:
		    *   0  1  2  3
		    *   4  5  6  7
		    *   8  9 10 11
		    *  12 13 14 15
		    */

public:

    /*! \brief Constructor for identity transformation.
     */
    Transformation() {
	_x[0] = _x[5] = _x[10] = _x[15] = 1.0; 
        _x[1] = _x[2] = _x[3] = _x[4] = _x[6] = _x[7] 
	    = _x[8] = _x[9] = _x[11] = _x[12] = _x[13]
	    = _x[14] = 0.0;
    }

    /*! \brief Constructor for preset transformation matrix.
     */
    Transformation( double x11, double x12, double x13, double x14,
		    double x21, double x22, double x23, double x24,
		    double x31, double x32, double x33, double x34,
		    double x41, double x42, double x43, double x44 ) { 
        _x[0]  = x11;
        _x[1]  = x12;
        _x[2]  = x13;
        _x[3]  = x14;
        _x[4]  = x21;
        _x[5]  = x22;
        _x[6]  = x23;
        _x[7]  = x24;
        _x[8]  = x31;
        _x[9]  = x32;
        _x[10] = x33;
        _x[11] = x34;
        _x[12] = x41;
        _x[13] = x42;
        _x[14] = x43;
        _x[15] = x44;
    }

    /*! \brief Copy constructor.
     */
    Transformation( const Transformation &m ) { 
        memcpy( _x, m._x, 16*sizeof(double) );
    }

    /*! \brief Destructor.
     */
    ~Transformation() {}





    /*! \brief Indexing for transformation matrix.
     */
    double &operator[]( int i ) {
        return( _x[i] );
    }

    /*! \brief Indexing for constant transformation matrix.
     */
    const double &operator[]( int i ) const {
        return( _x[i] );
    }
    
    /*! \brief Return determinant of matrix.
     */
    double determinant( void ) const;

    /*! \brief Return inverse matrix.
     */
    Transformation inverse( void ) const;

    /*! \brief Multiplication of transformation matrices for combining
     *  transformations.
     *
     *  Tranformation is done by multiplying the matrix with a vector
     *  from the right. Therefore the multiplication of transformation
     *  matrices has the effect that the right-hand-side
     *  transformation is applied first and left-hand-side second.
     */
    Transformation operator*( const Transformation &m ) const;





    /*! \brief Transform point \a xin.
     */
    Vec3D transform_point( const Vec3D &xin ) const;

    /*! \brief Inverse transform point \a xin.
     *
     *  This is a convenience function to inverting a transformation
     *  matrix and then doing a transform. If more than one transform
     *  is done inverse() and transform() functions should be used.
     */
    Vec3D inv_transform_point( const Vec3D &xin ) const;




    /*! \brief Transform vector \a xin.
     */
    Vec3D transform_vector( const Vec3D &xin ) const;

    /*! \brief Inverse transform vector \a xin.
     *
     *  This is a convenience function to inverting a transformation
     *  matrix and then doing a transform. If more than one transform
     *  is done inverse() and transform_vector() functions should be used.
     */
    Vec3D inv_transform_vector( const Vec3D &xin ) const;




    /*! \brief Translate transformation.
     *
     *  The effect of the new transformation is to first translate the
     *  coordinates and then apply the old transformation.
     */
    void translate( const Vec3D &d ) {
	Transformation t1 = translation( d );
	Transformation t2 = *this * t1;
	*this = t2;
    }

    /*! \brief Scale transformation.
     *
     *  The effect of the new transformation is to first scale the
     *  coordinates and then apply the old transformation.
     */
    void scale( const Vec3D &s ) {
	Transformation t1 = scaling( s );
	Transformation t2 = *this * t1;
	*this = t2;
    }

    /*! \brief Rotate transformation around x-axis.
     *
     *  The effect of the new transformation is to first rotate the
     *  coordinates and then apply the old transformation.
     */
    void rotate_x( double a ) {
	Transformation t1 = rotation_x( a );
	Transformation t2 = *this * t1;
	*this = t2;
    }

    /*! \brief Rotate transformation around y-axis.
     *
     *  The effect of the new transformation is to first rotate the
     *  coordinates and then apply the old transformation.
     */
    void rotate_y( double a ) {
	Transformation t1 = rotation_y( a );
	Transformation t2 = *this * t1;
	*this = t2;
    }

    /*! \brief Rotate transformation around z-axis.
     *
     *  The effect of the new transformation is to first rotate the
     *  coordinates and then apply the old transformation.
     */
    void rotate_z( double a ) {
	Transformation t1 = rotation_z( a );
	Transformation t2 = *this * t1;
	*this = t2;
    }



    /*! \brief Return translation transformation.
     */
    static Transformation translation( const Vec3D &d ) {
	return( Transformation(  1,  0,  0, d[0],
				 0,  1,  0, d[1],
				 0,  0,  1, d[2],
				 0,  0,  0,    1 ) );
    }

    /*! \brief Return scaling transformation.
     */
    static Transformation scaling( const Vec3D &s ) {
	return( Transformation( s[0],    0,    0,  0,
				   0, s[1],    0,  0,
				   0,    0, s[2],  0,
				   0,    0,    0,  1 ) );
    }
 
    /*! \brief Return rotation transformation rotating around x-axis.
     */
    static Transformation rotation_x( double a ) {
	return( Transformation(  1,      0,       0,  0,
				 0, cos(a), -sin(a),  0,
				 0, sin(a),  cos(a),  0,
				 0,      0,       0,  1 ) );
    }

    /*! \brief Return rotation transformation rotating around y-axis.
     */
    static Transformation rotation_y( double a ) {
	return( Transformation(  cos(a),  0, sin(a),  0,
			              0,  1,      0,  0,
			        -sin(a),  0, cos(a),  0,
			              0,  0,      0,  1 ) );
    }

    /*! \brief Return rotation transformation rotating around z-axis.
     */
    static Transformation rotation_z( double a ) {
	return( Transformation( cos(a), -sin(a),  0,  0,
				sin(a),  cos(a),  0,  0,
				     0,       0,  1,  0,
				     0,       0,  0,  1 ) );
    }


};




#endif

