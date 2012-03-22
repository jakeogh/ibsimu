/*! \file stl_solid.cpp
 *  \brief Stereolithography CAD file handling
 */

/* Copyright (c) 2011-2012 Taneli Kalvas. All rights reserved.
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
#include <ctype.h>
#include <string.h>
#include "stlfile.hpp"


//#define DEBUG_STL 1
#define EPS 1.0e-9


bool ciscomp( const char *str1, const char *str2, size_t n )
{
    size_t str1len = strlen( str1 );
    size_t str2len = strlen( str2 );

    if( str1len < n ) {
	if( str1len != str2len )
	    return( false );
	n = str1len;
    }
    if( str2len < n )
	return( false );

    for( size_t i = 0; i < n ; i++ ) {
        if( tolower(str1[i]) != tolower(str2[i]) )
            return( false );
    }
    return( true );
}


void STLFile::Triangle::read_binary_float_vector( Vec3D &x, std::ifstream &ifstr )
{
    for( size_t a = 0; a < 3; a++ ) {
	float r;
	ifstr.read( (char *)&r, 4 );
	x[a] = r;
    }
}


void STLFile::Triangle::read_ascii_float_vector( Vec3D &x, const char *buf, const std::string &filename, int linec )
{
    char *endptr;

    // Skip preceding whitespace
    while( isspace( *buf ) )
	buf++;

    for( size_t i = 0; i < 3; i++ ) {
	x[i] = strtod( buf, &endptr );
	if( endptr == buf ) {
	    throw( Error( ERROR_LOCATION, "unexpected end of line reading file \'" 
			  + filename + "\' on line " + to_string(linec) ) );
	}
	buf = endptr;
	while( isspace( *buf ) )
	    buf++;
    }

}


STLFile::Triangle::Triangle( std::ifstream &ifstr )
{
    read_binary_float_vector( _normal, ifstr );
    read_binary_float_vector( _p[0], ifstr );
    read_binary_float_vector( _p[1], ifstr );
    read_binary_float_vector( _p[2], ifstr );

    uint16_t attr;
    ifstr.read( (char *)&attr, 2 );
    _attr = attr;
}


STLFile::Triangle::Triangle( std::ifstream &ifstr, const char *buf, const std::string &filename, int &linec )
{
    //std::cout << "  Making new Triangle\n";

    if( !ciscomp( buf, "facet normal", 12 ) )
	throw( Error( ERROR_LOCATION, "Unexpected input on line " + to_string(linec) +
		      ", expecting \'facet normal\'." ) );
    buf += 12;
    read_ascii_float_vector( _normal, buf, filename, linec );

    // Tag "outer loop"
    std::string str;
    std::getline( ifstr, str );
    buf = str.c_str();
    linec++;

    // Skip whitespace
    while( isspace( *buf ) )
	buf++;

    if( !ciscomp( buf, "outer loop", 10 ) )
	throw( Error( ERROR_LOCATION, "Unexpected input on line " + to_string(linec) +
		      ", expecting \'outer loop\'." ) );

    // Read 3 vertices
    for( size_t i = 0; i < 3; i++ ) {
	std::getline( ifstr, str );
	buf = str.c_str();
	linec++;
	
	// Skip whitespace
	while( isspace( *buf ) )
	    buf++;
	
	if( !ciscomp( buf, "vertex", 6 ) )
	    throw( Error( ERROR_LOCATION, "Unexpected input on line " + to_string(linec) +
			  ", expecting \'vertex\'." ) );

	buf += 6;
	read_ascii_float_vector( _p[i], buf, filename, linec );
    }

    // Tag "endloop"
    std::getline( ifstr, str );
    buf = str.c_str();
    linec++;

    // Skip whitespace
    while( isspace( *buf ) )
	buf++;

    if( !ciscomp( buf, "endloop", 7 ) )
	throw( Error( ERROR_LOCATION, "Unexpected input on line " + to_string(linec) +
		      ", expecting \'endloop\'." ) );

    // Tag "endfacet"
    std::getline( ifstr, str );
    buf = str.c_str();
    linec++;

    // Skip whitespace
    while( isspace( *buf ) )
	buf++;

    if( !ciscomp( buf, "endfacet", 8 ) )
	throw( Error( ERROR_LOCATION, "Unexpected input on line " + to_string(linec) +
		      ", expecting \'endfacet\'." ) );


}


STLFile::Triangle::~Triangle()
{

}


const Vec3D &STLFile::Triangle::operator[]( int i ) const
{
    if( i == 0 )
	return( _p[0] );
    else if( i == 1 )
	return( _p[1] );
    else if( i == 2 )
	return( _p[2] );
    else
	throw( Error( ERROR_LOCATION, "Indexing error" ) );
}


int STLFile::Triangle::ray_cross( const Vec3D &x, const Vec3D &dir ) const
{
#ifdef DEBUG_STL
    std::cout << "ray_cross( " << x[0] << ", " << x[1] << ", " << x[2] << " )\n";
#endif

    // Find vectors for two edges sharing _p[0] and normal
    Vec3D u = _p[1] - _p[0];
    Vec3D v = _p[2] - _p[0];
    Vec3D n = cross( u, v );
    if( n == 0.0 )
	throw( Error( ERROR_LOCATION, "Zero triangle area" ) );	

    Vec3D w0 = x-_p[0];
    double a = -n*w0;
    double b = n*dir;
    if( fabs(b) < 1.0e-6 ) { // Ray parallel to triangle plane
	if( a == 0.0 ) // Ray in triangle plane
	    return( 2 );
	else
	    return( 0 ); // No collision
    }

    double r = a/b;
    if( r < 0.0 ) // Ray goes away from triangle
	return( 0 ); // No collision

    Vec3D I = x + r*dir; // Intersection point
    
    double uu = u*u;
    double uv = u*v;
    double vv = v*v;
    Vec3D w = I - _p[0];
    double wu = w*u;
    double wv = w*v;
    double D = uv*uv - uu*vv;

    double s = (uv * wv - vv * wu) / D;
    if( s < 0.0-EPS || s > 1.0+EPS ) 
	return( 0 );
    else if( s < 0.0+EPS || s > 1.0-EPS ) 
	return( 2 ); // Not sure, do perturbation
    double t = (uv * wu - uu * wv ) / D;
    if( t < 0.0-EPS || (s+t) > 1.0+EPS )
	return( 0 );
    else if( t < 0.0+EPS || (s+t) > 1.0-EPS )
	return( 2 ); // Not sure, do perturbation

    return( 1 ); // Collision
}


const Vec3D &STLFile::Triangle::normal( void ) const
{
    return( _normal );
}


const Vec3D &STLFile::Triangle::p1( void ) const
{
    return( _p[0] );
}


const Vec3D &STLFile::Triangle::p2( void ) const
{
    return( _p[1] );
}


const Vec3D &STLFile::Triangle::p3( void ) const
{
    return( _p[2] );
}




void STLFile::Triangle::debug_print( std::ostream &os ) const
{
    os << "**Triangle\n";
    os << "  normal = " << _normal << "\n";
    os << "  p1     = " << _p[0] << "\n";
    os << "  p2     = " << _p[1] << "\n";
    os << "  p3     = " << _p[2] << "\n";
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
    bbox_ppoint( min, max, _p[0] );
    bbox_ppoint( min, max, _p[1] );
    bbox_ppoint( min, max, _p[2] );
}



/* ******************** *
 * VTriangle            *
 * ******************** */


STLFile::VTriangle::VTriangle( uint32_t v1, uint32_t v2, uint32_t v3, const Vec3D &normal )
{
    _v[0] = v1;
    _v[1] = v2;
    _v[2] = v3;
    _normal = normal;
}


STLFile::VTriangle::VTriangle( const uint32_t v[3], const Vec3D &normal )
{
    _v[0] = v[0];
    _v[1] = v[1];
    _v[2] = v[2];
    _normal = normal;
}


STLFile::VTriangle::~VTriangle()
{
    
}


const Vec3D &STLFile::VTriangle::normal( void ) const
{
    return( _normal );
}


void STLFile::VTriangle::debug_print( std::ostream &os ) const
{
    os << "**VTriangle\n";    
    os << "  v = "  
       << std::setw(6) << _v[0] << " "
       << std::setw(6) << _v[1] << " "
       << std::setw(6) << _v[2] << "\n";
    os << "  normal = " << _normal << "\n";
}


const uint32_t &STLFile::VTriangle::operator[]( int i ) const
{
    return( _v[i] );
}


std::ostream &operator<<( std::ostream &os, const STLFile::VTriangle &vtri ) 
{
    os << std::setw(6) << to_string(vtri[0]) << " ";
    os << std::setw(6) << to_string(vtri[1]) << " ";
    os << std::setw(6) << to_string(vtri[2]) << " ";
    return( os );
}



/* ******************** *
 * STLFile              *
 * ******************** */


void STLFile::read_binary( std::ifstream &ifstr )
{
    // Read number of triangles
    uint32_t tcount;
    ifstr.read( (char *)(&tcount), 4 );

    // Check if sensible
    if( tcount > 10000000 )
	throw( ErrorUnimplemented( ERROR_LOCATION, "Too high number of triangles" ) );
    _triangle.reserve( tcount );
    
    for( uint32_t a = 0; a < tcount; a++ ) {
	_triangle.push_back( Triangle( ifstr ) );
    }
}


void STLFile::read_ascii( std::ifstream &ifstr )
{
    //std::cout << "Reading ASCII STL\n";

    // Read first line (header)
    std::string header;
    std::getline( ifstr, header );

    //std::cout << "  header = \'" + header + "\'\n";

    // Read line-by-line
    int linec = 0;
    while( !ifstr.eof() ) {
        
        std::string str;
        std::getline( ifstr, str );
	//std::cout << "  line = \'" + str + "\'\n";
	const char *buf = str.c_str();
        linec++;

	// Skip whitespace
	while( isspace( *buf ) )
	    buf++;

	if( ciscomp( buf, "endsolid", 8 ) )
	    break;

	_triangle.push_back( Triangle( ifstr, buf, _filename, linec ) );
    }

    //std::cout << "  Done\n";

    //throw( ErrorUnimplemented( ERROR_LOCATION, "Ascii STL file reader unimplemented" ) );
}


STLFile::STLFile( const std::string &filename )
{
    _filename = filename;
    std::ifstream ifstr( filename.c_str(), std::ios_base::binary );
    if( !ifstr.good() )
	throw( Error( ERROR_LOCATION, "Couldn't open file \'" + filename + "\'" ) );

    // Check if ascii
    // Ascii files start with "solid "
    // Binary files have free form 80 byte header
    char buf[80];
    ifstr.read( buf, 80 );
    if( !strncasecmp( buf, "solid ", 6 ) ) {
	_ascii = true;
	ifstr.seekg( 0 );
	read_ascii( ifstr );
    } else {
	_ascii = false;
	read_binary( ifstr );
    }

    ifstr.close();

    //build_vtriangle_data();
}


STLFile::~STLFile()
{
    
}


void STLFile::build_vtriangle_data( void )
{
    _vtri.clear();
    _vertex.clear();
    _vtri.reserve( _triangle.size() );

    // Go through all triangles and add coordinates to vertex list if
    // they are not already there. Add vertex triangles in the same
    // time using made vertices.
    for( size_t a = 0; a < _triangle.size(); a++ ) {

	uint32_t vert[3]; // Vertex indexes for triangle

	// Go through three vertices of triangle
	for( size_t vi = 0; vi < 3; vi++ ) {

	    // Search for vertex indexes, add vertices of triangle if
	    // they don's exist already.
	    Vec3D v = _triangle[a][vi];
	    uint32_t b;
	    for( b = 0; b < _vertex.size(); b++ ) {
		if( norm2(v-_vertex[b]) < EPS )
		    // Vertex found
		    break;
	    }
	    if( b == _vertex.size() ) {
		// New vertex needed
		_vertex.push_back( v );
	    }
	    vert[vi] = b;
	}

	// Add triangle using vertices
	_vtri.push_back( VTriangle( vert, _triangle[a].normal() ) );
    }
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


/*
// Tester using vtri and vertex arrays and normals
bool STLFile::inside( const Vec3D &x, double eps ) const
{
    //std::cout << "inside()\n";

    // Search for closest vertex to x
    int32_t imin = -1;
    double dmin = std::numeric_limits<double>::infinity();
    for( uint32_t a = 0; a < _vertex.size(); a++ ) {
	
	double d = norm2( x - _vertex[a] );
	if( d < dmin ) {
	    dmin = d;
	    imin = a;
	}
    }

    //std::cout << "x         = " << x << "\n";
    //std::cout << "imin      = " << imin << "\n";
    //std::cout << "dmin      = " << dmin << "\n";
    //std::cout << "vertex    = " << _vertex[imin] << "\n";

    // No vertices found
    if( imin == -1 )
	return( false );

    // Build a list of triangles connected to vertex imin
    std::vector<uint32_t> tri;
    //std::cout << "triangles = ";
    for( size_t a = 0; a < _vtri.size(); a++ ) {

	for( size_t b = 0; b < 3; b++ ) {
	    if( _vtri[a][b] == (uint32_t)imin ) {
		//std::cout << a << " ";
		tri.push_back( a );
		break;
	    }
	}
    }
    //std::cout << "\n";

    // Check if inside of all triangles on the list.
    // Return false if point is outside even one triangle.
    for( size_t a = 0; a < tri.size(); a++ ) {

	//std::cout << "test " << a << ": triangle " << tri[a] << ":\n";
	uint32_t vertex = _vtri[tri[a]][0];
	const Vec3D &p1 = _vertex[vertex];
	const Vec3D &n  = _vtri[tri[a]].normal();
	double pxn = (x-p1)*n;
	//std::cout << "  p1    = " << p1 << "\n";
	//std::cout << "  x-p1  = " << x-p1 << "\n";
	//std::cout << "  n     = " << n << "\n";
	//std::cout << "  pxn   = " << pxn << "\n";
	if( pxn > 0.0 ) {
	    //std::cout << "OUT\n";
	    return( false );
	}
    }

    //std::cout << "IN\n";
    return( true );
}
*/


bool STLFile::inside( const Vec3D &x, double eps ) const
{
#ifdef DEBUG_STL
    std::cout << "\n\ninside( " << x[0] << ", " << x[1] << ", " << x[2] << " )\n";
#endif

    Vec3D dir( 1.0, 0.0, 0.0 );
    Vec3D X( x );

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
	    stat = t->ray_cross( X, dir );
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
	    X[1] += 1.0e-6;
	else
	    X[2] += 1.0e-6;
    }

    std::stringstream se;
    se << "Perturbation failed\n"
       << "  x = " << x << "\n"
       << "  triangle no = " << a << "\n"
       << "  tri.p1 = " << t->p1() << "\n"
       << "  tri.p2 = " << t->p2() << "\n"
       << "  tri.p3 = " << t->p3() << "\n";
    throw( Error( ERROR_LOCATION, se.str() ) );
}


void STLFile::debug_print( std::ostream &os ) const
{
    os << "**STLFile\n";
    os << "  trianglec = " << _triangle.size() << "\n";
    
    for( size_t a = 0; a < _triangle.size(); a++ ) {
	os << "  triangle[" << a << "]:\n";
	_triangle[a].debug_print( os );
    }

    os << "  vertexc = " << _vertex.size() << "\n";
    for( size_t a = 0; a < _vertex.size(); a++ )
	os << "  vertex[" << a << "] = " << _vertex[a] << "\n";

    os << "  vtric = " << _vtri.size() << "\n";
    for( size_t a = 0; a < _vtri.size(); a++ ) {
	os << "  vtri[" << a << "] = " << _vtri[a] << ", ";
	os << "normal = " << _vtri[a].normal() << "\n";
    }
}


