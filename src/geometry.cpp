/*! \file geometry.cpp
 *  \brief Source code for geometry.cpp
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

#include <iostream>
#include <iomanip>
#include <math.h>
#include "geometry.hpp"
#include "func_solid.hpp"
#include "error.hpp"
#include "verbose.hpp"
#include "file.hpp"


Geometry::Geometry( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h )
    : _geom_mode(geom_mode), _size(size), _origo(origo)
{
    _h = fabs(h);
    if( _h == 0.0 )
	throw( Error( ERROR_LOCATION, "zero mesh step size" ) );

    if( verbose_output ) {
	Int3D one(1,1,1);
	std::cout << "Constructing geometry\n";
	std::cout << "  origo = " << origo << "\n";
	std::cout << "  size  = " << size << "\n";
	std::cout << "  max   = " << origo+_h*(size-one) << "\n";
	std::cout << "  h     = " << _h << "\n";
    }

    if( _geom_mode == MODE_CYL ) {
	if( _origo[1] < 0.0 )
	    throw( Error( ERROR_LOCATION, "negative origo in r-direction" ) );
    }
    
    if( _geom_mode == MODE_3D ) {
	if( _size[0] < 3 || _size[1] < 3 || _size[2] < 3 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    } else if( _geom_mode == MODE_2D || _geom_mode == MODE_CYL ) {
	if( _size[0] < 3 || _size[1] < 3 || _size[2] != 1 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    } else {
	if( _size[0] < 3 || _size[1] != 1 || _size[2] != 1 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    }

    // Calculate vector max
    _max = Vec3D( _origo(0)+_h*(_size[0]-1),
		  _origo(1)+_h*(_size[1]-1),
		  _origo(2)+_h*(_size[2]-1) );

    _n = 0;
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );

    _built = false;
    _smesh = new signed char[_size[0]*_size[1]*_size[2]];
}


Geometry::Geometry( std::istream &s )
{
    _geom_mode = (geom_mode_e)read_int32( s );
    _size      = Int3D( s );
    _origo     = Vec3D( s );
    _h         = read_double( s );
    _n         = read_int32( s );

    // Calculate vector max
    _max = Vec3D( _origo(0)+_h*(_size[0]-1),
		  _origo(1)+_h*(_size[1]-1),
		  _origo(2)+_h*(_size[2]-1) );

    for( int32_t a = 0; a < _n; a++ ) {
	int32_t fileid = read_int32( s );
	if( fileid == FILEID_FUNCSOLID )
	    _sdata.push_back( new FuncSolid( s ) );
	else
	    throw( Error( ERROR_LOCATION, "unknown solid type" ) );
    }

    for( int32_t a = 0; a < _n+6; a++ )
	_bound.push_back( Bound( s ) );

    _built = read_int8( s );
    _smesh = new signed char[_size[0]*_size[1]*_size[2]];
    read_compressed_block( s, _size[0]*_size[1]*_size[2], (int8_t *)_smesh );
}


Geometry::~Geometry()
{
    for( int32_t a = 0; a < _n; a++ )
	delete _sdata[a];
    delete [] _smesh;
}


int32_t Geometry::dim( void ) const
{
    switch( _geom_mode ) {
    case MODE_1D:
	return( 1 );
	break;
    case MODE_2D:
	return( 2 );
	break;
    case MODE_CYL:
	return( 2 );
	break;
    default:
	return( 3 );
	break;
    }
}


void Geometry::set_solid( int32_t n, const Solid *s )
{
    if( n <= 6 || n > _n+7 )
	throw( Error( ERROR_LOCATION, "illegal solid number" + to_string(n) ) );

    if( n <= _n+6 ) {
	delete _sdata[n-7];
    } else {
	_sdata.push_back( 0 );
	_bound.push_back( Bound(BOUND_DIRICHLET,0.0) );
	_n++;
    }

    _sdata[n-7] = s;
}


const Solid *Geometry::get_solid( int32_t n ) const
{
    if( n <= 6 || n > _n+6 )
	throw( Error( ERROR_LOCATION, "illegal solid number " + to_string(n) ) );
    
    return( _sdata[n-7] );
}


void Geometry::set_boundary( int32_t n, const Bound &b )
{
    if( n <= 0 || n > _n+6 )
	throw( Error( ERROR_LOCATION, "illegal solid number " + to_string(n) ) );

    _bound[n-1] = b;
}


Bound Geometry::get_boundary( int32_t n ) const
{
    if( n <= 0 || n > _n+6 )
	throw( Error( ERROR_LOCATION, "illegal solid number " + to_string(n) ) );

    return( _bound[n-1] );
}


int32_t Geometry::inside( const Vec3D &x ) const
{
    for( ssize_t a = _sdata.size()-1; a >= 0 ; a-- ) {
	if( _sdata[a]->inside( x ) )
	    return( a+7 );
    }
    if( x[2] > _max[2] )
	return( 6 );
    else if( x[2] < _origo[2] )
	return( 5 );
    else if( x[1] > _max[1] )
	return( 4 );
    else if( x[1] < _origo[1] )
	return( 3 );
    else if( x[0] > _max[0] )
	return( 2 );
    else if( x[0] < _origo[0] )
	return( 1 );

    return( 0 );
}


bool Geometry::inside( int32_t n, const Vec3D &x ) const
{
    if( n <= 6 || n > _n+6 )
	throw( Error( ERROR_LOCATION, "illegal solid number n=" + to_string(n) ) );

    return( _sdata[n-7]->inside( x ) );
}


signed char Geometry::mesh_check( int32_t i, int32_t j, int32_t k ) const
{
    if( i < 0 )
	return( 1 );
    else if( i >= _size[0] )
	return( 2 );
    if( j < 0 )
	return( 3 );
    else if( j >= _size[1] )
	return( 4 );
    if( k < 0 )
	return( 5 );
    else if( k >= _size[2] )
	return( 6 );

    return( _smesh[i + j*_size[0] + k*_size[0]*_size[1]] );
}


double Geometry::bracket_surface( int32_t n, const Vec3D &xin, const Vec3D &xout, Vec3D &xsurf ) const
{
    int a;
    Vec3D xl = xin;
    Vec3D xh = xout;

    // Do iteration
    for( a = 0; a < 12; a++ ) {
	xsurf = 0.5*(xl+xh);
	if( inside( n, xsurf ) )
	    xl = xsurf;
	else
	    xh = xsurf;
    }

    // Calculate best guess, the midpoint
    xsurf = 0.5*(xl+xh);

    // Return parametric distance
    for( a = 0; a < 3; a++ ) {
	if( xin[a] != xout[a] )
	    return( (xsurf[a] - xin[a]) / (xout[a] - xin[a]) );
    }

    // If endpoints are the same, return zero
    return( 0.0 );
}


void Geometry::build_mesh( void )
{
    int a;
    int32_t i, j, k;
    double x, y, z;

    _built = true;

    // Set solids
    for( k = 0; k < _size[2]; k++ ) {
	z = k*_h+_origo[2];
	for( j = 0; j < _size[1]; j++ ) {
	    y = j*_h+_origo[1];
	    for( i = 0; i < _size[0]; i++ ) {
		x = i*_h+_origo[0];
		mesh(i,j,k) = inside( Vec3D(x,y,z) );
	    }
	}
    }

    // Mark solid edges
    for( k = 0; k < _size[2]; k++ ) {
	for( j = 0; j < _size[1]; j++ ) {
	    for( i = 0; i < _size[0]; i++ ) {
		// Check if node is an edge node
		if( (a = mesh(i,j,k)) != 0 &&
		    ( mesh_check(i-1,j,  k  ) == 0 ||
		      mesh_check(i+1,j,  k  ) == 0 ||
		      mesh_check(i,  j-1,k  ) == 0 || 
		      mesh_check(i,  j+1,k  ) == 0 ||
		      mesh_check(i,  j,  k-1) == 0 || 
		      mesh_check(i,  j,  k+1) == 0 ) ) {
		    // Mark node as an edge
		    mesh(i,j,k) = -a;
		}
	    }
	}
    }

    // Mark boundaries
    // Mark xmin
    if( _bound[0].type == BOUND_NEUMANN ) a = -1;
    else a = 1;
    for( k = 0; k < _size[2]; k++ ) {
	for( j = 0; j < _size[1]; j++ ) {
	    if( mesh(0,j,k) == 0 )
		mesh(0,j,k) = a;
	}
    }
    // Mark xmax
    if( _bound[1].type == BOUND_NEUMANN ) a = -2;
    else a = 2;
    for( k = 0; k < _size[2]; k++ ) {
	for( j = 0; j < _size[1]; j++ ) {
	    if( mesh(_size[0]-1,j,k) == 0 )
		mesh(_size[0]-1,j,k) = a;
	}
    }
    if( _geom_mode == MODE_2D || _geom_mode == MODE_CYL ||
	_geom_mode == MODE_3D ) {
	// Mark ymin.
	if( _bound[2].type == BOUND_NEUMANN ) a = -3;
	else a = 3;
	for( k = 0; k < _size[2]; k++ ) {
	    for( i = 0; i < _size[0]; i++ ) {
		if( mesh(i,0,k) == 0 )
		    mesh(i,0,k) = a;
	    }
	}
	// Mark ymax
	if( _bound[3].type == BOUND_NEUMANN ) a = -4;
	else a = 4;
	for( k = 0; k < _size[2]; k++ ) {
	    for( i = 0; i < _size[0]; i++ ) {
		if( mesh(i,_size[1]-1,k) == 0 )
		    mesh(i,_size[1]-1,k) = a;
	    }
	}
    }
    if( _geom_mode == MODE_3D ) {
	// Mark zmin
	if( _bound[4].type == BOUND_NEUMANN ) a = -5;
	else a = 5;
	for( j = 0; j < _size[1]; j++ ) {
	    for( i = 0; i < _size[0]; i++ ) {
		if( mesh(i,j,0) == 0 )
		    mesh(i,j,0) = a;
	    }
	}
	// Mark ymax
	if( _bound[5].type == BOUND_NEUMANN ) a = -6;
	else a = 6;
	for( j = 0; j < _size[1]; j++ ) {
	    for( i = 0; i < _size[0]; i++ ) {
		if( mesh(i,j,_size[2]-1) == 0 )
		    mesh(i,j,_size[2]-1) = a;
	    }
	}
    }
}


void Geometry::save( std::ostream &s ) const
{
    write_int32( s, _geom_mode );
    _size.save( s );
    _origo.save( s );
    write_double( s, _h );
    write_int32( s, _n );
    for( int32_t a = 0; a < _n; a++ )
	_sdata[a]->save( s );
    for( int32_t a = 0; a < _n+6; a++ )
	_bound[a].save( s );

    write_int8( s, _built );
    write_compressed_block( s, _size[0]*_size[1]*_size[2]*sizeof(signed char), _smesh );
}


void Geometry::debug_print( void ) const
{
    std::cout << "**Geometry\n";

    if( _geom_mode == MODE_1D )
	std::cout << "geom_mode = MODE_1D\n";
    else if( _geom_mode == MODE_2D )
	std::cout << "geom_mode = MODE_2D\n";
    else if( _geom_mode == MODE_CYL )
	std::cout << "geom_mode = MODE_CYL\n";
    else if( _geom_mode == MODE_3D )
	std::cout << "geom_mode = MODE_3D\n";
    else
	std::cout << "geom_mode = Unknown\n";
    std::cout << "size = (" 
	      << _size[0] << ", "
	      << _size[1] << ", "
	      << _size[2] << ")\n";
    std::cout << "origo = (" 
	      << _origo[0] << ", "
	      << _origo[1] << ", "
	      << _origo[2] << ")\n";
    std::cout << "max = (" 
	      << _max[0] << ", "
	      << _max[1] << ", "
	      << _max[2] << ")\n";
    std::cout << "h = " << _h << "\n";
    std::cout << "n = " << _n << "\n";
    if( _n == 0 )
	std::cout << "no sdata\n";
    for( int32_t a = 0; a < _n; a++ ) {
	std::cout << "sdata[" << a << "]:\n";
	_sdata[a]->debug_print();
    }
    for( int32_t a = 0; a < _n+6; a++ ) {
	std::cout << "bound[" << a << "] = (";
	if( _bound[a].type == BOUND_NEUMANN )
	    std::cout << "BOUND_NEUMANN, ";
	else 
	    std::cout << "BOUND_DIRICHLET, ";
	std::cout << _bound[a].val << ")\n";
    }
    std::cout << "built = " << _built << "\n";
    std::cout << "smesh = \n";
    for( int32_t j = 0; j < _size[1]; j++ ) {
	for( int32_t i = 0; i < _size[0]; i++ )
	    std::cout << std::setw(2) << (int)mesh(i,j,0) << " ";
	std::cout << "\n";
    }
}













