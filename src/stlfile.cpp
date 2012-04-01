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
#include "ibsimu.hpp"


#define DEBUG_STL 1


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


/* ******************** *
 * Triangle             *
 * ******************** */


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


uint16_t STLFile::Triangle::attr( void ) const
{
    return( _attr );
}


const Vec3D &STLFile::Triangle::normal( void ) const
{
    return( _normal );
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


STLFile::VTriangle::VTriangle( uint32_t v1, uint32_t v2, uint32_t v3 )
{
    _v[0] = v1;
    _v[1] = v2;
    _v[2] = v3;
}


STLFile::VTriangle::VTriangle( const uint32_t v[3] )
{
    _v[0] = v[0];
    _v[1] = v[1];
    _v[2] = v[2];
}


STLFile::VTriangle::~VTriangle()
{
    
}


void STLFile::VTriangle::debug_print( std::ostream &os ) const
{
    os << "**VTriangle\n";    
    os << "  v = "  
       << std::setw(6) << _v[0] << " "
       << std::setw(6) << _v[1] << " "
       << std::setw(6) << _v[2] << "\n";
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
    _tri.reserve( tcount );
    
    for( uint32_t a = 0; a < tcount; a++ ) {
	_tri.push_back( Triangle( ifstr ) );
    }
}


void STLFile::read_ascii( std::ifstream &ifstr )
{
    // Read first line (header)
    std::string header;
    std::getline( ifstr, header );

    // Read line-by-line
    int linec = 0;
    while( !ifstr.eof() ) {
        
        std::string str;
        std::getline( ifstr, str );
	const char *buf = str.c_str();
        linec++;

	// Skip whitespace
	while( isspace( *buf ) )
	    buf++;

	if( ciscomp( buf, "endsolid", 8 ) )
	    break;

	_tri.push_back( Triangle( ifstr, buf, _filename, linec ) );
    }
}


STLFile::STLFile( const std::string &filename,
		  double vertex_matching_eps, 
		  double signed_volume_eps )
    : _filename(filename), _vertex_matching_eps(vertex_matching_eps), 
      _signed_volume_eps(signed_volume_eps)
{
    ibsimu.message( 1 ) << "Reading STL-file \'" << filename << "\'\n";
    ibsimu.inc_indent();

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

    build_vtriangle_data();
    check_vtriangle_data();
    calculate_bbox();
    offset_vtriangle_data();
    _tri.clear();
    _vpos.assign( _vertex.size(), false );
    _vneg.assign( _vertex.size(), false );

    ibsimu.dec_indent();
}


STLFile::~STLFile()
{
    
}


void STLFile::offset_vtriangle_data( void )
{
    _offset = -_bbox[0] + 1.0e-3*(_bbox[1]-_bbox[0]);

    for( size_t a = 0; a < _vertex.size(); a++ )
	_vertex[a] += _offset;
}


void STLFile::build_vtriangle_data( void )
{
    ibsimu.message(1) << "Making vertex connections\n";

    _vtri.clear();
    _vertex.clear();
    _vtri.reserve( _tri.size() );

    // Go through all triangles and add coordinates to vertex list if
    // they are not already there. Add vertex triangles in the same
    // time using made vertices.

    for( size_t a = 0; a < _tri.size(); a++ ) {

	uint32_t vert[3]; // Vertex indexes for triangle

	// Go through three vertices of triangle
	for( size_t vi = 0; vi < 3; vi++ ) {

	    // Search for vertex indexes, add vertices of triangle if
	    // they don's exist already.
	    Vec3D v = _tri[a][vi];
	    uint32_t b;
	    for( b = 0; b < _vertex.size(); b++ ) {
		if( norm2(v-_vertex[b]) < _vertex_matching_eps )
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
	_vtri.push_back( VTriangle( vert ) );
    }
}


void STLFile::check_vtriangle_data( void )
{
    ibsimu.message( 1 ) << "Checking mesh data\n";
    ibsimu.inc_indent();

    ibsimu.message( 1 ) << _vtri.size() << " triangles\n";
    ibsimu.message( 1 ) << _vertex.size() << " vertices\n";

    // Check 
    // 1. Triangle edge pairing. Every triangle must have exactly 
    // one neighbouring triangle for each of the three edges.
    // 2. Neighbouring triangles must be defined in same direction.
    // 3. Zero area triangles.
    for( uint32_t a = 0; a < _vtri.size(); a++ ) {

	Vec3D V1 = _vertex[_vtri[a][1]] - _vertex[_vtri[a][0]];
	Vec3D V2 = _vertex[_vtri[a][2]] - _vertex[_vtri[a][0]];
	double area = 0.5*norm2( cross( V1, V2 ) );
	if( fabs(area) == 0.0 )
	    throw( Error( ERROR_LOCATION, "Zero area triangle " + to_string(a) ) );
	for( uint32_t b = 0; b < 3; b++ ) {
	    
	    int e1 = _vtri[a][b];
	    int e2 = _vtri[a][(b+1)%3];

	    bool found = false;
	    uint32_t neighbour;
	    uint32_t c;
	    for( c = 0; c < _vtri.size(); c++ ) {
		if( c == a ) continue; // No self-checking
		for( uint32_t d = 0; d < 3; d++ ) {
		    
		    int f1 = _vtri[c][d];
		    int f2 = _vtri[c][(d+1)%3];
		    if( f1 == e1 && f2 == e2 ) {
			// Incorrect orientation
			throw( Error( ERROR_LOCATION, "Incorrect orientation between neighbouring triangles " +
				      to_string(a) + " and " + to_string(c) ) );
		    } else if( f1 == e2 && f2 == e1 ) {
			// Correct orientation
			if( found == true )
			    throw( Error( ERROR_LOCATION, "Double neighbours (" + to_string(neighbour) +
					  " and " + to_string(c) + ") found for triangle " +
					  to_string(a) ) );
			found = true;
			neighbour = c;
		    }
		}
	    }
	    if( !found ) {
		throw( Error( ERROR_LOCATION, "Triangle " + to_string(a) + " neighbour not found" ) );
	    }
	}
    }

    ibsimu.dec_indent();    
}


size_t STLFile::vertexc( void )
{
    return( _vertex.size() );
}

size_t STLFile::trianglec( void )
{
    return( _vtri.size() );
}


void STLFile::get_bbox( Vec3D &min, Vec3D &max ) const
{
    min = _bbox[0];
    max = _bbox[1];
}


void STLFile::calculate_bbox( void )
{
    Vec3D min = Vec3D( std::numeric_limits<double>::infinity(),
		       std::numeric_limits<double>::infinity(),
		       std::numeric_limits<double>::infinity() );
    Vec3D max = Vec3D( -std::numeric_limits<double>::infinity(),
		       -std::numeric_limits<double>::infinity(),
		       -std::numeric_limits<double>::infinity() );

    for( uint32_t a = 0; a < _tri.size(); a++ ) {
	_tri[a].update_bbox( min, max );
    }

    _bbox[0] = min;
    _bbox[1] = max;
}


#define STL_EDGE1   0
#define STL_EDGE2   1
#define STL_EDGE3   2
#define STL_FACE    3
#define STL_OUTSIDE 4
#define STL_INSIDE  5


int STLFile::signvol4( const Vec3D &q0, const Vec3D &q1, 
		       const Vec3D &q2, const Vec3D &q3 )
{
    double x = q0[0]*(-q1[1]*q2[2] + q1[1]*q3[2] + q2[1]*q1[2] - 
		       q2[1]*q3[2] - q3[1]*q1[2] + q3[1]*q2[2] ) +
	       q1[0]*( q0[1]*q2[2] - q0[1]*q3[2] - q2[1]*q0[2] + 
		       q2[1]*q3[2] + q3[1]*q0[2] - q3[1]*q2[2] ) +
	       q2[0]*(-q0[1]*q1[2] + q0[1]*q3[2] + q1[1]*q0[2] - 
		       q1[1]*q3[2] - q3[1]*q0[2] + q3[1]*q1[2] ) +
	       q3[0]*( q0[1]*q1[2] - q0[1]*q2[2] - q1[1]*q0[2] + 
		       q1[1]*q2[2] + q2[1]*q0[2] - q2[1]*q1[2] );
    if( x < _signed_volume_eps ) {
	if( x <= -_signed_volume_eps )
	    return( -1 );
	return( 0 );
    }
    return( 1 );
}


int STLFile::signvol3( const Vec3D &q1, const Vec3D &q2, const Vec3D &q3 )
{
    double x = q1[0]*( q2[1]*q3[2] - q3[1]*q2[2] ) -
	       q1[1]*( q2[0]*q3[2] - q3[0]*q2[2] ) +
	       q1[2]*( q2[0]*q3[1] - q3[0]*q2[1] );
    if( x < _signed_volume_eps ) {
	if( x <= -_signed_volume_eps )
	    return( -1 );
	return( 0 );
    }
    return( 1 );
}


/* Classify inclusion of point p in tetrahedron (o,q1,q2,q3), when
 * The sense of tetrahedron (o,q1,q2,q3) is ss.
 */
int STLFile::classify_original_tetrahedron( int ss, const Vec3D &p, 
					    const Vec3D &q1, const Vec3D &q2, const Vec3D &q3 )
{
    int s0,s1,s2,s3;

    if( ss > 0 ) {
	// Positive sense (o,q1,q2,q3)
	if( (s1 = signvol3( p, q2, q3 )) < 0 )
	    return( STL_OUTSIDE );
	if( (s2 = signvol3( p, q3, q1 )) < 0 )
	    return( STL_OUTSIDE );
	if( (s3 = signvol3( p, q1, q2 )) < 0 )
	    return( STL_OUTSIDE );
	if( (s0 = signvol4( p, q1, q2, q3 )) < 0 )
	    return( STL_OUTSIDE );
    } else {
	// Negative sense (o,q1,q2,q3)
	if( (s1 = signvol3( p, q2, q3 )) > 0 )
	    return( STL_OUTSIDE );
	if( (s2 = signvol3( p, q3, q1 )) > 0 )
	    return( STL_OUTSIDE );
	if( (s3 = signvol3( p, q1, q2 )) > 0 )
	    return( STL_OUTSIDE );
	if( (s0 = signvol4( p, q1, q2, q3 )) > 0 )
	    return( STL_OUTSIDE );
    }

    // Edges and faces
    if( s1 == 0 ) {
	if( s2 == 0 )
	    return( STL_EDGE3 );
	else if( s3 == 0 )
	    return( STL_EDGE2 );
	return( STL_FACE );
    } else if( s2 == 0 ) {
	if( s3 == 0 )
	    return( STL_EDGE1 );
	return( STL_FACE );
    } else if( s3 == 0 ) {
	return( STL_FACE );
    }

    return( STL_INSIDE );
}


bool STLFile::inside( const Vec3D &p )
{
    Vec3D x = p+_offset;
    for( uint32_t a = 0; a < 3; a++ ) {
	if( x[a] <= 0.0 )
	    return( false );
    }

    // Clear positive and negative vertex arrays
    _vpos.assign( _vertex.size(), false );
    _vneg.assign( _vertex.size(), false );

    int incl = 0;
    for( uint32_t a = 0; a < _vtri.size(); a++ ) {

	int ss = signvol3( _vertex[_vtri[a][0]],
			   _vertex[_vtri[a][1]],
			   _vertex[_vtri[a][2]] );
	int stat = classify_original_tetrahedron( ss, x,
						  _vertex[_vtri[a][0]],
						  _vertex[_vtri[a][1]],
						  _vertex[_vtri[a][2]] );
	if( stat == STL_INSIDE ) {
	    incl += 2*ss;
	} else if( stat == STL_FACE ) {
	    incl += ss;
	} else if( stat != STL_OUTSIDE ) {
	    if( ss > 0 && !_vpos[_vtri[a][stat]] ) {
		_vpos[_vtri[a][stat]] = true;
		incl += 2*ss;
	    } else if( ss < 0 && !_vneg[_vtri[a][stat]] ) {
		_vneg[_vtri[a][stat]] = true;
		incl += 2*ss;
	    }
	}
    }

    if( incl > 0 )
	return( true ); 
    return( false );
}


void STLFile::debug_print( std::ostream &os ) const
{
    os << "**STLFile\n";

    os << "  vertexc = " << _vertex.size() << "\n";
    for( size_t a = 0; a < _vertex.size(); a++ )
	os << "  vertex[" << a << "] = " << _vertex[a] << "\n";

    os << "  vtric = " << _vtri.size() << "\n";
    for( size_t a = 0; a < _vtri.size(); a++ ) {
	os << "  vtriangle[" << a << "]:\n";
	os << "    v1 = " << _vtri[a][0] << "\n";
	os << "    v2 = " << _vtri[a][1] << "\n";
	os << "    v3 = " << _vtri[a][2] << "\n";
    }
}


