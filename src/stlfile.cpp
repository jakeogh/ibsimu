/*! \file stl_solid.cpp
 *  \brief Stereolithography CAD file handling
 */

/* Copyright (c) 2011 Taneli Kalvas. All rights reserved.
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

#include <cstddef>
#include <limits>
#include <sstream>
#include "stlfile.hpp"


//#define DEBUG_STL 1


void STLFile::Triangle::read_binary_float_vector( Vec3D &x, std::ifstream &ifstr )
{
    for( size_t a = 0; a < 3; a++ ) {
	float r;
	ifstr.read( (char *)&r, 4 );
	x[a] = r;
    }
}


STLFile::Triangle::Triangle( std::ifstream &ifstr )
{
    read_binary_float_vector( _normal, ifstr );
    read_binary_float_vector( _p1, ifstr );
    read_binary_float_vector( _p2, ifstr );
    read_binary_float_vector( _p3, ifstr );

    uint16_t attr;
    ifstr.read( (char *)&attr, 2 );
    _attr = attr;
}


STLFile::Triangle::~Triangle()
{

}


int STLFile::Triangle::ray_cross( const Vec3D &x, const Vec3D &l ) const
{
#ifdef DEBUG_STL
    std::cout << "ray_cross( " << x[0] << ", " << x[1] << ", " << x[2] << " )\n";
#endif

    // Find vectors for two edges sharing _p1
    Vec3D edge1 = _p2 - _p1;
    Vec3D edge2 = _p3 - _p1;

    // begin calculating determinant - also used to calculate U parameter
    Vec3D pvec = cross( l, edge2 );
    //Vec3D pvec( 0.0, -edge2[2], edge2[1] );

    // if determinant is near zero, ray lies in plane of triangle
    double det = edge1 * pvec;

    // calculate distance from p1 to ray origin x
    Vec3D tvec = x - _p1;
    double inv_det = 1.0 / det;

    Vec3D qvec;
#ifdef DEBUG_STL
    std::cout << "det = " << det << "\n";
#endif
    if( det > 1.0e-6 ) {

        // calculate U parameter and test bounds
        double u = tvec * pvec;
#ifdef DEBUG_STL
	std::cout << "u = " << u << "\n";
#endif
        if( u < 0.0 || u > det ) {
#ifdef DEBUG_STL
	    std::cout << "no intersection\n";
#endif
            return( 0 );
	}

        // prepare to test V parameter
        qvec = cross( tvec, edge1 );

        // calculate V parameter and test bounds
        double v = l*qvec;
#ifdef DEBUG_STL
	std::cout << "v = " << v << "\n";
#endif
        if( v < 0.0 || u + v > det ) {
#ifdef DEBUG_STL
	    std::cout << "no intersection\n";
#endif
            return( 0 );
	}

    } else if( det < -1.0e-6 ) {

        // calculate U parameter and test bounds
        double u = tvec * pvec;
#ifdef DEBUG_STL
	std::cout << "u = " << u << "\n";
#endif
        if( u > 0.0 || u < det ) {
#ifdef DEBUG_STL
	    std::cout << "no intersection\n";
#endif
            return( 0 );
	}

        // prepare to test V parameter
        qvec = cross( tvec, edge1 );

        // calculate V parameter and test bounds
        double v = l*qvec;
#ifdef DEBUG_STL
	std::cout << "v = " << v << "\n";
#endif
        if( v > 0.0 || u + v < det ) {
#ifdef DEBUG_STL
	    std::cout << "no intersection\n";
#endif
            return( 0 );
	}

    } else {

        // ray is parallel to the plane of the triangle
#ifdef DEBUG_STL
	std::cout << "parallel and in plane, go for perturbation\n";
#endif
	return( 2 );

    }

    // ray intersects triangle, calculate a
    double a = (edge2 * qvec) * inv_det;
#ifdef DEBUG_STL
    std::cout << "a = " << a << "\n";
#endif
    if( a > 0 ) {
#ifdef DEBUG_STL
	std::cout << "no intersection\n";
#endif
	return( 0 );
    }

#ifdef DEBUG_STL
    std::cout << "intersection\n";
#endif
    return( 1 );
}


const Vec3D &STLFile::Triangle::normal( void ) const
{
    return( _normal );
}


const Vec3D &STLFile::Triangle::p1( void ) const
{
    return( _p1 );
}


const Vec3D &STLFile::Triangle::p2( void ) const
{
    return( _p2 );
}


const Vec3D &STLFile::Triangle::p3( void ) const
{
    return( _p3 );
}




void STLFile::Triangle::debug_print( std::ostream &os ) const
{
    os << "**Triangle\n";
    os << "  normal = " << _normal << "\n";
    os << "  p1     = " << _p1 << "\n";
    os << "  p2     = " << _p2 << "\n";
    os << "  p3     = " << _p3 << "\n";
    os << "  attr   = " << _attr << "\n";
}


void STLFile::Triangle::bbox_ppoint( Vec3D &min, Vec3D &max, const Vec3D &p )
{
    for( int a = 0; a < 3; a++ ) {
	if( p[a] < min[a] )
	    min[a] = p[a];
	if( p[a] > max[a] )
	    max[a] = p[a];
    }
}


void STLFile::Triangle::update_bbox( Vec3D &min, Vec3D &max ) const
{
    bbox_ppoint( min, max, _p1 );
    bbox_ppoint( min, max, _p2 );
    bbox_ppoint( min, max, _p3 );
}


void STLFile::read_binary( std::ifstream &ifstr )
{
    // Read number of triangles
    uint32_t tcount;
    ifstr.read( (char *)(&tcount), 4 );

    _triangle.reserve( tcount );
    
    for( uint32_t a = 0; a < tcount; a++ ) {
	_triangle.push_back( Triangle( ifstr ) );
    }
}


STLFile::STLFile( const std::string &filename )
{
    std::ifstream ifstr( filename.c_str() );
    if( !ifstr.good() )
	throw Error( ERROR_LOCATION, "Couldn't open file \'" + filename + "\'" );

    // Check if ascii
    // Ascii files start with "solid "
    // Binary files have free form 80 byte header
    char buf[80];
    ifstr.read( buf, 80 );
    if( !strncmp( buf, "Solid ", 6 ) ) {
	_ascii = true;
	ifstr.seekg( 0 );
	throw( ErrorUnimplemented( ERROR_LOCATION, "Ascii STL file reader unimplemented" ) );
    } else {
	_ascii = false;
	read_binary( ifstr );
    }

    ifstr.close();
}


STLFile::~STLFile()
{
    
}


size_t STLFile::size( void )
{
    return( _triangle.size() );
}


void STLFile::get_bbox( Vec3D &min, Vec3D &max ) const
{
    min = Vec3D( std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity() );
    max = Vec3D( -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity(),
		 -std::numeric_limits<double>::infinity() );

    for( uint32_t a = 0; a < _triangle.size(); a++ ) {
	_triangle[a].update_bbox( min, max );
    }
}


bool STLFile::inside( const Vec3D &x, double eps ) const
{
#ifdef DEBUG_STL
    std::cout << "\n\ninside( " << x[0] << ", " << x[1] << ", " << x[2] << " )\n";
#endif

    Vec3D l( 1.0, 0.0, 0.0 );

    // Perturbation loop 
    uint32_t a;
    const Triangle *t = NULL;
    for( uint32_t b = 0; b < 3; b++ ) {

	int stat = 0;
        int par = 0;
        for( a = 0; a < _triangle.size(); a++ ) {

#ifdef DEBUG_STL
	    std::cout << "triangle " << a << ":\n";
#endif

	    t = &_triangle[a];
	    stat = t->ray_cross( x, l );
	    if( stat == 1 )
		par = !par;
	    else if( stat == 2 ) 
		break;
        }

        if( stat != 2 )
            return( par );

	// Perturbation
#ifdef DEBUG_STL
	std::cout << "do perturbation\n";
#endif
	if( b == 0 )
	    l[1] += 0.216;
	else
	    l[2] += 0.117;
    }

    std::stringstream se;
    se << "Perturbation failed\n"
       << "  x = " << x << "\n"
       << "  triangle no = " << a << "\n"
       << "  tri.p1 = " << t->p1() << "\n"
       << "  tri.p2 = " << t->p2() << "\n"
       << "  tri.p3 = " << t->p3() << "\n";
    throw Error( ERROR_LOCATION, se.str() );
}


void STLFile::debug_print( std::ostream &os ) const
{
    os << "**STLFile\n";
    os << "  trianglec = " << _triangle.size() << "\n";
    
    for( size_t a = 0; a < _triangle.size(); a++ )
	_triangle[a].debug_print( os );
}
