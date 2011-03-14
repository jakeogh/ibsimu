/*! \file geometry.cpp
 *  \brief %Geometry definition.
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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

#include <iostream>
#include <iomanip>
#include <math.h>
#include <string.h>
#include "geometry.hpp"
#include "func_solid.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "file.hpp"


std::ostream &operator<<( std::ostream &os, const Bound &b )
{
    if( b.type == BOUND_NEUMANN )
	os << "(BOUND_NEUMANN, ";
    else 
	os << "(BOUND_DIRICHLET, ";
    os << b.val << ")";
    return( os );
}


Geometry::Geometry( geom_mode_e geom_mode, Int3D size, Vec3D origo, double h )
    : Mesh(geom_mode,size,origo,h)
{
    check_definition();

    if( ibsimu.get_verbose_output() ) {
	std::cout << "Constructing geometry\n";
	std::cout << "  origo     = " << _origo << "\n";
	std::cout << "  size      = " << _size << "\n";
	std::cout << "  max       = " << _max << "\n";
	std::cout << "  h         = " << _h << "\n";
	std::cout << "  nodecount = " << nodecount() << "\n";
    }

    _n = 0;
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );
    _bound.push_back( Bound(BOUND_NEUMANN,0.0) );

    _built = false;
    _smesh = new uint32_t[_size[0]*_size[1]*_size[2]];
}


Geometry::Geometry( std::istream &s )
    : Mesh(s)
{
    check_definition();

    _n = read_int32( s );
    for( uint32_t a = 0; a < _n; a++ ) {
	int32_t fileid = read_int32( s );
	if( fileid == FILEID_FUNCSOLID )
	    _sdata.push_back( new FuncSolid( s ) );
	else
	    throw( Error( ERROR_LOCATION, "unknown solid type" ) );
    }

    for( uint32_t a = 0; a < _n+6; a++ )
	_bound.push_back( Bound( s ) );

    _built = read_int8( s );
    _smesh = new uint32_t[_size[0]*_size[1]*_size[2]];
    read_compressed_block( s, sizeof(uint32_t)*_size[0]*_size[1]*_size[2], 
			   (int8_t *)_smesh );

    uint32_t nearsolidsize = read_int32( s );
    _nearsolid.resize( nearsolidsize );
    read_compressed_block( s, sizeof(uint8_t)*nearsolidsize, 
			   (int8_t *)&_nearsolid[0] );
}


Geometry::~Geometry()
{
    for( uint32_t a = 0; a < _n; a++ )
	delete _sdata[a];
    delete [] _smesh;
}


void Geometry::check_definition()
{
    // Geometry specific checks
    if( _geom_mode == MODE_3D ) {
	if( _size[0] < 3 || _size[1] < 3 || _size[2] < 3 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    } else if( _geom_mode == MODE_2D || _geom_mode == MODE_CYL ) {
	if( _geom_mode == MODE_CYL ) {
	    if( _origo(1) != 0.0 )
		throw( Error( ERROR_LOCATION, "j=0 node not on axis" ) );
	}
	if( _size[0] < 3 || _size[1] < 3 || _size[2] != 1 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    } else {
	if( _size[0] < 3 || _size[1] != 1 || _size[2] != 1 )
	    throw( Error( ERROR_LOCATION, "illegal mesh size" ) );
    }
}


void Geometry::set_solid( uint32_t n, const Solid *s )
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


uint32_t Geometry::number_of_solids() const
{
    return( _n );
}


uint32_t Geometry::number_of_boundaries() const
{
    return( _n+6 );
}


const Solid *Geometry::get_solid( uint32_t n ) const
{
    if( n <= 6 || n > _n+6 )
	throw( Error( ERROR_LOCATION, "illegal solid number " + to_string(n) ) );
    
    return( _sdata[n-7] );
}


void Geometry::set_boundary( uint32_t n, const Bound &b )
{
    if( n <= 0 || n > _n+6 )
	throw( Error( ERROR_LOCATION, "illegal solid number " + to_string(n) ) );

    if( n >= 7 && b.type != BOUND_DIRICHLET )
	throw( Error( ERROR_LOCATION, "trying to set solid " + to_string(n) + " as Neumann boundary" ) );

    _bound[n-1] = b;
}


Bound Geometry::get_boundary( uint32_t n ) const
{
    if( n <= 0 || n > _n+6 )
	throw( Error( ERROR_LOCATION, "illegal solid number " + to_string(n) ) );

    return( _bound[n-1] );
}


uint32_t Geometry::inside( const Vec3D &x ) const
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


bool Geometry::inside( uint32_t n, const Vec3D &x ) const
{
    if( n <= 6 || n > _n+6 )
	throw( Error( ERROR_LOCATION, "illegal solid number n=" + to_string(n) ) );

    return( _sdata[n-7]->inside( x ) );
}


uint32_t Geometry::mesh_check( int32_t i, int32_t j, int32_t k ) const
{
    if( i < 0 )
	return( SMESH_NODE_ID_DIRICHLET | 1 );
    else if( i >= _size[0] )
	return( SMESH_NODE_ID_DIRICHLET | 2 );
    if( j < 0 )
	return( SMESH_NODE_ID_DIRICHLET | 3 );
    else if( j >= _size[1] )
	return( SMESH_NODE_ID_DIRICHLET | 4 );
    if( k < 0 )
	return( SMESH_NODE_ID_DIRICHLET | 5 );
    else if( k >= _size[2] )
	return( SMESH_NODE_ID_DIRICHLET | 6 );

    return( _smesh[i + j*_size[0] + k*_size[0]*_size[1]] );
}


double Geometry::bracket_surface( uint32_t n, const Vec3D &xin, const Vec3D &xout, Vec3D &xsurf ) const
{
    Vec3D xl = xin;
    Vec3D xh = xout;

    // Do iteration
    for( uint32_t a = 0; a < 8; a++ ) {
	xsurf = 0.5*(xl+xh);
	if( inside( n, xsurf ) )
	    xl = xsurf;
	else
	    xh = xsurf;
    }

    // Calculate best guess, the midpoint
    xsurf = 0.5*(xl+xh);

    // Return parametric distance
    for( uint32_t a = 0; a < 3; a++ ) {
	if( xin[a] != xout[a] )
	    return( (xsurf[a] - xin[a]) / (xout[a] - xin[a]) );
    }

    // If endpoints are the same, return zero
    return( 0.0 );
}


uint8_t Geometry::bracket_ndist( int32_t i, int32_t j, int32_t k, int32_t solid, int sign, int coord ) const
{
    // Floating point version of bracketing
    /*
    Vec3D vout( _origo[0]+i*_h, _origo[1]+j*_h, _origo[2]+k*_h );
    double xsurf = 0.5;
    double xin = 1.0;
    double xout = 0.0;

    // Do iteration
    for( uint32_t a = 0; a < 8; a++ ) {
	Vec3D vtest( vout );
	vtest[coord] += sign*xsurf*_h;
	if( inside( solid, vtest ) )
	    xin = xsurf;
	else
	    xout = xsurf;
	xsurf = 0.5*(xin+xout);
    }
    */

    Vec3D vout( _origo[0]+i*_h, _origo[1]+j*_h, _origo[2]+k*_h );
    uint32_t bp = 0x80;
    uint32_t surf = 0x80;

    // Do iteration
    double step = sign*_h;
    for( uint32_t a = 0; a < 8; a++ ) {
	Vec3D vtest( vout );
	vtest[coord] += step*surf/255.0;
	if( inside( solid, vtest ) )
	    surf -= bp;
	bp = bp >> 1;
	surf += bp;
    }

    return( surf );
}


uint32_t Geometry::is_solid( int32_t i, int32_t j, int32_t k ) const
{
    uint32_t a = mesh_check( i, j, k );
    if( (a & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET &&
	(a & SMESH_BOUNDARY_NUMBER_MASK) >= 7 )
	return( a & SMESH_BOUNDARY_NUMBER_MASK );
    return( 0 );
}


bool Geometry::is_near_solid( int32_t i, int32_t j, int32_t k ) const
{
    return( is_solid(i-1,j,  k  ) ||
	    is_solid(i+1,j,  k  ) ||
	    is_solid(i,  j-1,k  ) || 
	    is_solid(i,  j+1,k  ) ||
	    is_solid(i,  j,  k-1) || 
	    is_solid(i,  j,  k+1) );
}


void Geometry::add_near_solid_distance( std::vector<uint8_t> &ndist, uint8_t dist )
{
    // zero distance is illegal, can cause problems with solver
    if( dist == 0 )
	ndist.push_back( 1 );
    else
	ndist.push_back( dist );
}


void Geometry::add_near_solid_entry( uint32_t &near_solid_index, int32_t i, int32_t j, int32_t k )
{
    uint32_t solid = 0;         // Solid number of neighbour
    int nsolids = 0;            // Number of near neighbour solids nodes.
    uint8_t neighbours = 0;     // Bit flags for neighbours
    std::vector<uint8_t> ndist; // Neighbour distances

    // X
    if( (solid = is_solid(i-1,j,k)) ) {
	neighbours += 0x20;
	nsolids++;
	add_near_solid_distance( ndist, bracket_ndist( i,j,k, solid, -1, 0 ) );
    }
    neighbours = neighbours >> 1;
    if( ( solid = is_solid(i+1,j,k)) ) {
	neighbours += 0x20;
	nsolids++;
	add_near_solid_distance( ndist, bracket_ndist( i,j,k, solid, +1, 0 ) );
    }
    neighbours = neighbours >> 1;

    // Y
    if( (solid = is_solid(i,j-1,k)) ) {
	neighbours += 0x20;
	nsolids++;
	add_near_solid_distance( ndist, bracket_ndist( i,j,k, solid, -1, 1 ) );
    }
    neighbours = neighbours >> 1;
    if( (solid = is_solid(i,j+1,k)) ) {
	neighbours += 0x20;
	nsolids++;
	add_near_solid_distance( ndist, bracket_ndist( i,j,k, solid, +1, 1 ) );
    }
    neighbours = neighbours >> 1;

    // Z
    if( (solid = is_solid(i,j,k-1)) ) {
	neighbours += 0x20;
	nsolids++;
	add_near_solid_distance( ndist, bracket_ndist( i,j,k, solid, -1, 2 ) );
    }
    neighbours = neighbours >> 1;
    if( (solid = is_solid(i,j,k+1)) ) {
	neighbours += 0x20;
	nsolids++;
	add_near_solid_distance( ndist, bracket_ndist( i,j,k, solid, +1, 2 ) );
    }

    size_t ind = _nearsolid.size();
    near_solid_index = ind + 1 + sizeof(uint8_t)*nsolids;
    _nearsolid.resize( near_solid_index );
    _nearsolid[ind] = neighbours;
    memcpy( (void *)&_nearsolid[ind+1], (void *)&ndist[0], sizeof(uint8_t)*nsolids );
}


void Geometry::build_mesh( void )
{
    if( ibsimu.get_verbose_output() )
	std::cout << "Building mesh\n";

    _built = true;

    // Mark solid (Dirichlet) nodes. Others left to zero.
    for( int32_t k = 0; k < _size[2]; k++ ) {
	double z = k*_h+_origo[2];
	for( int32_t j = 0; j < _size[1]; j++ ) {
	    double y = j*_h+_origo[1];
	    for( int32_t i = 0; i < _size[0]; i++ ) {
		double x = i*_h+_origo[0];
		uint32_t nid = inside( Vec3D(x,y,z) );
		if( nid )
		    mesh(i,j,k) = SMESH_NODE_ID_DIRICHLET | nid;
		else
		    mesh(i,j,k) = 0;
	    }
	}
    }

    // Mark 0 nodes on boundaries as Dirichlet, Neumann boundaries or
    // as near solid nodes. Here x overrides y, which overrides
    // z. Order is important as near solid indexes are processed also
    // by the solver preprocess/postprocess routines.
    uint32_t nid;
    uint32_t near_solid_index = 0;

    // Mark xmin and xmax
    for( uint32_t bound = 1; bound <= 2; bound++ ) {
	int32_t i = 0;
	if( bound == 2 ) i = _size[0]-1;
	if( get_boundary(bound).type == BOUND_NEUMANN )
	    nid = SMESH_NODE_ID_NEUMANN | bound;
	else
	    nid = SMESH_NODE_ID_DIRICHLET | bound;
	for( int32_t k = 0; k < _size[2]; k++ ) {
	    for( int32_t j = 0; j < _size[1]; j++ ) {
		if( mesh(i,j,k) == 0 ) {
		    if( get_boundary(bound).type == BOUND_NEUMANN && is_near_solid(i,j,k) ) {
			// Mark node as near solid vacuum and build near solid data
			mesh(i,j,k) = SMESH_NODE_ID_NEAR_SOLID | near_solid_index;
			add_near_solid_entry( near_solid_index, i, j, k );
		    } else {
			mesh(i,j,k) = nid;
		    }
		}
	    }
	}
    }

    if( _geom_mode == MODE_2D || _geom_mode == MODE_CYL ||
	_geom_mode == MODE_3D ) {
	// Ymin and Ymax
	for( uint32_t bound = 3; bound <= 4; bound++ ) {
	    int32_t j = 0;
	    if( bound == 4 ) j = _size[1]-1;
	    if( get_boundary(bound).type == BOUND_NEUMANN )
		nid = SMESH_NODE_ID_NEUMANN | bound;
	    else
		nid = SMESH_NODE_ID_DIRICHLET | bound;
	    for( int32_t k = 0; k < _size[2]; k++ ) {
		for( int32_t i = 0; i < _size[0]; i++ ) {
		    if( mesh(i,j,k) == 0 ) {
			if( get_boundary(bound).type == BOUND_NEUMANN && is_near_solid(i,j,k) ) {
			    // Mark node as near solid vacuum and build near solid data
			    mesh(i,j,k) = SMESH_NODE_ID_NEAR_SOLID | near_solid_index;
			    add_near_solid_entry( near_solid_index, i, j, k );
			} else {
			    mesh(i,j,k) = nid;
			}
		    }
		}
	    }
	}
    }
    if( _geom_mode == MODE_3D ) {
	// Zmin and Zmax
	for( uint32_t bound = 5; bound <= 6; bound++ ) {
	    int32_t k = 0;
	    if( bound == 6 ) k = _size[2]-1;
	    if( get_boundary(bound).type == BOUND_NEUMANN )
		nid = SMESH_NODE_ID_NEUMANN | bound;
	    else
		nid = SMESH_NODE_ID_DIRICHLET | bound;
	    for( int32_t j = 0; j < _size[1]; j++ ) {
		for( int32_t i = 0; i < _size[0]; i++ ) {
		    if( mesh(i,j,k) == 0 ) {
			if( get_boundary(bound).type == BOUND_NEUMANN && is_near_solid(i,j,k) ) {
			    // Mark node as near solid vacuum and build near solid data
			    mesh(i,j,k) = SMESH_NODE_ID_NEAR_SOLID | near_solid_index;
			    add_near_solid_entry( near_solid_index, i, j, k );
			} else {
			    mesh(i,j,k) = nid;
			}
		    }
		}
	    }
	}
    }

    // Mark 0 nodes next to solid nodes as near solid nodes and build
    // nearsolid data array, rest of the 0 nodes are marked as pure
    // vacuum. No conflict should occur because 0 node = pure vacuum.
    for( int32_t k = 0; k < _size[2]; k++ ) {
	for( int32_t j = 0; j < _size[1]; j++ ) {
	    for( int32_t i = 0; i < _size[0]; i++ ) {
		if( mesh(i,j,k) != 0 )
		    continue;
		if( is_near_solid(i,j,k) ) {
		    // Mark node as near solid vacuum and build near solid data
		    mesh(i,j,k) = SMESH_NODE_ID_NEAR_SOLID | near_solid_index;
		    add_near_solid_entry( near_solid_index, i, j, k );
		} else {
		    // Pure vacuum
		    mesh(i,j,k) = SMESH_NODE_ID_PURE_VACUUM;
		}
	    }
	}
    }

    // Report node counts
    if( ibsimu.get_verbose_output() ) {
	int b;
	uint32_t ncount = nodecount();
	uint32_t nvacuum = 0;
	uint32_t nnearsolid = 0;
	uint32_t nneumann = 0;
	uint32_t ndirichlet = 0;
	int nsolid[_n];
	for( uint32_t a = 0; a < _n; a++ )
	    nsolid[a] = 0;

	for( uint32_t a = 0; a < ncount; a++ ) {

	    switch( mesh(a) & SMESH_NODE_ID_MASK ) {
	    case SMESH_NODE_ID_NEAR_SOLID:
		nnearsolid++;
		break;
	    case SMESH_NODE_ID_PURE_VACUUM:
		nvacuum++;
		break;
	    case SMESH_NODE_ID_NEUMANN:
		nneumann++;
		if( (b = (mesh(a) & SMESH_BOUNDARY_NUMBER_MASK)) >= 7 )
		    nsolid[b-7]++;
		break;
	    case SMESH_NODE_ID_DIRICHLET:
		ndirichlet++;
		if( (b = (mesh(a) & SMESH_BOUNDARY_NUMBER_MASK)) >= 7 )
		    nsolid[b-7]++;
		break;
	    }
	}
	
	std::cout << "  Done. Built mesh with:\n";
	std::cout << "  " << nvacuum << " pure vacuum nodes\n";
	std::cout << "  " << nnearsolid << " near solid nodes\n";
	std::cout << "  " << nneumann << " neumann nodes\n";
	std::cout << "  " << ndirichlet << " dirichlet nodes\n";
	for( uint32_t a = 0; a < _n; a++ )
	    std::cout << "  " << nsolid[a] << " solid " << a+7 << " nodes\n";
    }
}


uint8_t Geometry::solid_dist( uint32_t i, uint32_t j, uint32_t k, uint32_t dir ) const
{
    uint32_t snode = _smesh[i + j*_size[0] + k*_size[0]*_size[1]];
    if( (snode & SMESH_NODE_ID_MASK) != SMESH_NODE_ID_NEAR_SOLID )
	throw( Error( ERROR_LOCATION, "not a near solid node" ) );

    const uint8_t *nptr = &_nearsolid[snode & SMESH_NEAR_SOLID_INDEX_MASK];
    uint8_t neighbours = nptr[0];
    nptr++;
    uint32_t a = 0;
    while( a < dir ) {
	if( neighbours & 0x01 )
	    nptr++;
	neighbours = neighbours >> 1;
	a++;
    }
    if( (neighbours & 0x01) == 0x00 )
	throw( Error( ERROR_LOCATION, (const std::string)"no near neighbour in selected direction"
		      ", dir = " + to_string(dir) 
		      + ", neighbours = " + to_string((int)_nearsolid[snode & SMESH_NEAR_SOLID_INDEX_MASK]) ) );

    return( *nptr );
}


uint8_t Geometry::solid_dist( uint32_t i, uint32_t dir ) const
{
    uint32_t snode = _smesh[i];
    if( (snode & SMESH_NODE_ID_MASK) != SMESH_NODE_ID_NEAR_SOLID )
	throw( Error( ERROR_LOCATION, "not a near solid node" ) );

    const uint8_t *nptr = &_nearsolid[snode & SMESH_NEAR_SOLID_INDEX_MASK];
    uint8_t neighbours = nptr[0];
    nptr++;
    uint32_t a = 0;
    while( a < dir ) {
	if( neighbours & 0x01 )
	    nptr++;
	neighbours = neighbours >> 1;
	a++;
    }
    if( (neighbours & 0x01) == 0x00 )
	throw( Error( ERROR_LOCATION, (const std::string)"no near neighbour in selected direction"
		      ", dir = " + to_string(dir) 
		      + ", neighbours = " + to_string((int)_nearsolid[snode & SMESH_NEAR_SOLID_INDEX_MASK]) ) );

    return( *nptr );
}


void Geometry::save( std::ostream &s ) const
{
    Mesh::save( s );
    write_int32( s, _n );
    for( uint32_t a = 0; a < _n; a++ )
	_sdata[a]->save( s );
    for( uint32_t a = 0; a < _n+6; a++ )
	_bound[a].save( s );

    write_int8( s, _built );
    write_compressed_block( s, _size[0]*_size[1]*_size[2]*sizeof(uint32_t), 
			    (int8_t *)_smesh );

    write_int32( s, _nearsolid.size() );
    write_compressed_block( s, _nearsolid.size()*sizeof(uint8_t), 
			    (int8_t *)&_nearsolid[0] );
}


void Geometry::debug_print( std::ostream &os ) const
{
    Mesh::debug_print( os );

    os << "**Geometry\n";

    os << "n = " << _n << "\n";
    if( _n == 0 )
	os << "no sdata\n";
    for( uint32_t a = 0; a < _n; a++ ) {
	os << "sdata[" << a << "]:\n";
	_sdata[a]->debug_print( os );
    }
    for( uint32_t a = 0; a < _n+6; a++ ) {
	os << "bound[" << a+1 << "] = " << _bound[a] << "\n";
    }
    os << "built = " << _built << "\n";
    if( (_geom_mode == MODE_1D || _geom_mode == MODE_2D || _geom_mode == MODE_CYL) && 
	_size[0] <= 20 && _size[1] <= 20 ) {
	os << "mesh visualization:\n";
	for( int32_t j = _size[1]-1; j >= 0; j-- ) {
	    for( int32_t i = 0; i < _size[0]; i++ ) {
		uint32_t ind = mesh(i,j);
		if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_PURE_VACUUM ) {
		    os << " ";
		} else if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEUMANN ) {
		    os << "N";
		} else if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET ) {
		    if( (ind & SMESH_BOUNDARY_NUMBER_MASK) <= 6 )
			os << "D";
		    else
			os << "S";
		} else if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEAR_SOLID ) {
		    os << ".";
		}
	    }
	    os << "\n";
	}
    }

    for( int32_t k = 0; k < _size[2]; k++ ) {
	for( int32_t j = 0; j < _size[1]; j++ ) {
	    for( int32_t i = 0; i < _size[0]; i++ ) {

		uint32_t ind = mesh(i,j,k);
		os << "smesh(" << i << ", " << j << ", " << k << ") = "
		   << "0x" << std::hex << std::setfill('0')
		   << std::setw(8) << ind << " ("
		   << std::dec << std::setfill(' ');

		if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_PURE_VACUUM ) {
		    os << "pure vacuum)\n";
		} else if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEUMANN ) {
		    os << "neumann, solid " << (ind & SMESH_BOUNDARY_NUMBER_MASK) << ")\n";
		} else if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_DIRICHLET ) {
		    os << "dirichlet, solid " << (ind & SMESH_BOUNDARY_NUMBER_MASK) << ")\n";
		} else if( (ind & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEAR_SOLID ) {		    
		    os << "near solid, index " << (ind & SMESH_NEAR_SOLID_INDEX_MASK) << ")\n";
		    ind = (ind & SMESH_NEAR_SOLID_INDEX_MASK);
		    uint8_t sflag = _nearsolid[ind];
		    os << std::hex << std::setfill('0')
		       << "  near solid flags = 0x" << std::setw(2) << (int)sflag << "\n"
		       << std::dec << std::setfill(' ');
		    uint8_t mask = 1;
		    uint8_t *ptr = (uint8_t *)&_nearsolid[ind+1];
		    for( uint32_t a = 0; a < 6; a++ ) {
			if( mask & sflag ) {
			    os << "  ndist[" << a << "] = " << (int)*ptr << "\n";
			    ptr++;
			}
			mask = mask << 1;
		    }
		}
	    }
	}
    }
}


