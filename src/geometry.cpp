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
#include <fstream>
#include <math.h>
#include "geometry.hpp"
#include "func_solid.hpp"
#include "dxf_solid.hpp"
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
	ibsimu.vout() << "Constructing geometry\n";
	ibsimu.vout() << "  origo     = " << _origo << "\n";
	ibsimu.vout() << "  size      = " << _size << "\n";
	ibsimu.vout() << "  max       = " << _max << "\n";
	ibsimu.vout() << "  h         = " << _h << "\n";
	ibsimu.vout() << "  nodecount = " << nodecount() << "\n";
    }

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


Geometry::Geometry( std::istream &is )
    : Mesh(is)
{
    check_definition();

    if( ibsimu.get_verbose_output() )
	ibsimu.vout() << "Constructing Geometry from stream\n";

    _n = read_int32( is );
    for( uint32_t a = 0; a < _n; a++ ) {
	int32_t fileid = read_int32( is );
	if( fileid == FILEID_FUNCSOLID )
	    _sdata.push_back( new FuncSolid( is ) );
	else if( fileid == FILEID_DXFSOLID )
	    _sdata.push_back( new DXFSolid( is ) );
	else
	    throw( Error( ERROR_LOCATION, "unknown solid type" ) );
    }

    for( uint32_t a = 0; a < _n+6; a++ )
	_bound.push_back( Bound( is ) );
    
    _built = read_int8( is );
    _smesh = new signed char[_size[0]*_size[1]*_size[2]];
    read_compressed_block( is, _size[0]*_size[1]*_size[2], (int8_t *)_smesh );
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
    double eps = 1.0e-6*h();
    if( x[2] > _max[2]+eps )
	return( 6 );
    else if( x[2] < _origo[2]-eps )
	return( 5 );
    else if( x[1] > _max[1]+eps )
	return( 4 );
    else if( x[1] < _origo[1]-eps )
	return( 3 );
    else if( x[0] > _max[0]+eps )
	return( 2 );
    else if( x[0] < _origo[0]-eps )
	return( 1 );

    return( 0 );
}


bool Geometry::inside( uint32_t n, const Vec3D &x ) const
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


double Geometry::bracket_surface( uint32_t n, const Vec3D &xin, const Vec3D &xout, Vec3D &xsurf ) const
{
    if( xin == xout ) 
	throw( Error( ERROR_LOCATION, "xin and xout are the same point" ) );

    Vec3D xl = xin;
    Vec3D xh = xout;

    // Do iteration
    for( int a = 0; a < 8; a++ ) {
	xsurf = 0.5*(xl+xh);
	if( inside( n, xsurf ) )
	    xl = xsurf;
	else
	    xh = xsurf;
    }

    // Calculate best guess, the midpoint
    xsurf = 0.5*(xl+xh);

    // Return parametric distance calculated using axis where
    // coordinate difference is largest
    int a;
    Vec3D dif = xout - xin;
    dif.abs();
    if( dif[0] > dif[1] ) {
 	if( dif[0] > dif[2] )
	    a = 0;
	else
	    a = 2;
    } else {
	if( dif[1] > dif[2] )
	    a = 1;
	else
	    a = 2;
    }
    return( (xsurf[a] - xin[a]) / (xout[a] - xin[a]) );
}


Vec3D Geometry::surface_normal( const Vec3D &x ) const
{
    switch( geom_mode() ) {
    case MODE_2D:
    case MODE_CYL:
	return( surface_normal_2d( x ) );
    case MODE_3D:
	return( surface_normal_3d( x ) );
    default:
	throw( ErrorUnimplemented( ERROR_LOCATION ) );
    }
}


Vec3D Geometry::surface_normal_3d( const Vec3D &x ) const
{
    double dx = 0.1*_h;

    Vec3D px[8];
    uint32_t p[8];
    int surc;
    bool sur[12];

    do {
	// Build corner coordinates
	px[0] = x+Vec3D( -dx, -dx, -dx );
	px[1] = x+Vec3D( -dx, +dx, -dx );
	px[2] = x+Vec3D( +dx, +dx, -dx );
	px[3] = x+Vec3D( +dx, -dx, -dx );
	px[4] = x+Vec3D( -dx, -dx, +dx );
	px[5] = x+Vec3D( -dx, +dx, +dx );
	px[6] = x+Vec3D( +dx, +dx, +dx );
	px[7] = x+Vec3D( +dx, -dx, +dx );

	// Check inside() in corners
	p[0] = inside( px[0] );
	p[1] = inside( px[1] );
	p[2] = inside( px[2] );
	p[3] = inside( px[3] );
	p[4] = inside( px[4] );
	p[5] = inside( px[5] );
	p[6] = inside( px[6] );
	p[7] = inside( px[7] );

	// Check if at least three intersections
	surc = 0;
	if( p[0] != p[1] ) {
	    surc++;
	    sur[0] = true;
	} else {
	    sur[0] = false;
	}
	if( p[1] != p[2] ) {
	    surc++;
	    sur[1] = true;
	} else {
	    sur[1] = false;
	}
	if( p[2] != p[3] ) {
	    surc++;
	    sur[2] = true;
	} else {
	    sur[2] = false;
	}
	if( p[3] != p[0] ) {
	    surc++;
	    sur[3] = true;
	} else {
	    sur[3] = false;
	}
	//
	if( p[4] != p[5] ) {
	    surc++;
	    sur[4] = true;
	} else {
	    sur[4] = false;
	}
	if( p[5] != p[6] ) {
	    surc++;
	    sur[5] = true;
	} else {
	    sur[5] = false;
	}
	if( p[6] != p[7] ) {
	    surc++;
	    sur[6] = true;
	} else {
	    sur[6] = false;
	}
	if( p[7] != p[4] ) {
	    surc++;
	    sur[7] = true;
	} else {
	    sur[7] = false;
	}
	//
	if( p[0] != p[4] ) {
	    surc++;
	    sur[8] = true;
	} else {
	    sur[8] = false;
	}
	if( p[1] != p[5] ) {
	    surc++;
	    sur[9] = true;
	} else {
	    sur[9] = false;
	}
	if( p[2] != p[6] ) {
	    surc++;
	    sur[10] = true;
	} else {
	    sur[10] = false;
	}
	if( p[3] != p[7] ) {
	    surc++;
	    sur[11] = true;
	} else {
	    sur[11] = false;
	}

	dx *= 2.0;
	if( dx >= _h )
	    return( Vec3D(0,0,0) );

    } while( surc < 3 );

    // Bracket surface locations
    int b = 0;
    Vec3D xs[12];
    if( sur[0] ) {
	if( p[0] )
	    bracket_surface( p[0], px[0], px[1], xs[b++] );
	else
	    bracket_surface( p[1], px[1], px[0], xs[b++] );
    }
    if( sur[1] ) {
	if( p[1] )
	    bracket_surface( p[1], px[1], px[2], xs[b++] );
	else
	    bracket_surface( p[2], px[2], px[1], xs[b++] );
    }
    if( sur[2] ) {
	if( p[2] )
	    bracket_surface( p[2], px[2], px[3], xs[b++] );
	else
	    bracket_surface( p[3], px[3], px[2], xs[b++] );
    }
    if( sur[3] ) {
	if( p[3] )
	    bracket_surface( p[3], px[3], px[0], xs[b++] );
	else
	    bracket_surface( p[0], px[0], px[3], xs[b++] );
    }
    //
    if( sur[4] ) {
	if( p[4] )
	    bracket_surface( p[4], px[4], px[5], xs[b++] );
	else
	    bracket_surface( p[5], px[5], px[4], xs[b++] );
    }
    if( sur[5] ) {
	if( p[5] )
	    bracket_surface( p[5], px[5], px[6], xs[b++] );
	else
	    bracket_surface( p[6], px[6], px[5], xs[b++] );
    }
    if( sur[6] ) {
	if( p[6] )
	    bracket_surface( p[6], px[6], px[7], xs[b++] );
	else
	    bracket_surface( p[7], px[7], px[6], xs[b++] );
    }
    if( sur[7] ) {
	if( p[7] )
	    bracket_surface( p[7], px[7], px[4], xs[b++] );
	else
	    bracket_surface( p[4], px[4], px[7], xs[b++] );
    }
    //
    if( sur[8] ) {
	if( p[0] )
	    bracket_surface( p[0], px[0], px[4], xs[b++] );
	else
	    bracket_surface( p[4], px[4], px[0], xs[b++] );
    }
    if( sur[9] ) {
	if( p[1] )
	    bracket_surface( p[1], px[1], px[5], xs[b++] );
	else
	    bracket_surface( p[5], px[5], px[1], xs[b++] );
    }
    if( sur[10] ) {
	if( p[2] )
	    bracket_surface( p[2], px[2], px[6], xs[b++] );
	else
	    bracket_surface( p[6], px[6], px[2], xs[b++] );
    }
    if( sur[11] ) {
	if( p[3] )
	    bracket_surface( p[3], px[3], px[7], xs[b++] );
	else
	    bracket_surface( p[7], px[7], px[3], xs[b++] );
    }

    // Choose for three first intersections for normal
    // Three points furthest from each other would be better, but more hassle to search
    Vec3D v1 = xs[1]-xs[0];
    Vec3D v2 = xs[2]-xs[0];
    Vec3D n = cross( v1, v2 );
    n.normalize();

    // Check polarity
    if( inside( x+0.1*h()*n ) )
	n *= -1.0;

    return( n );
}


Vec3D Geometry::surface_normal_2d( const Vec3D &x ) const
{
    double dx = 0.1*_h;

    Vec3D px[4];
    uint32_t p[4];
    int surc;
    bool sur[4];

    do {
	// Build corner coordinates
	px[0] = x+Vec3D(-dx,-dx,0);
	px[1] = x+Vec3D(-dx,+dx,0);
	px[2] = x+Vec3D(+dx,+dx,0);
	px[3] = x+Vec3D(+dx,-dx,0);

	// Check inside() in corners
	p[0] = inside( px[0] );
	p[1] = inside( px[1] );
	p[2] = inside( px[2] );
	p[3] = inside( px[3] );

	// Check if two intersections
	surc = 0;
	if( p[0] != p[1] ) {
	    surc++;
	    sur[0] = true;
	} else {
	    sur[0] = false;
	}
	if( p[1] != p[2] ) {
	    surc++;
	    sur[1] = true;
	} else {
	    sur[1] = false;
	}
	if( p[2] != p[3] ) {
	    surc++;
	    sur[2] = true;
	} else {
	    sur[2] = false;
	}
	if( p[3] != p[0] ) {
	    surc++;
	    sur[3] = true;
	} else {
	    sur[3] = false;
	}

	dx *= 2.0;
	if( dx >= _h )
	    return( Vec3D(0,0,0) );

    } while( surc != 2 );

    // Bracket surface locations
    int pol = 0;
    int b = 0;
    Vec3D xs[2];
    for( int a = 0; a < 4; a++ ) {
	if( sur[a] ) {
	    pol = p[a];
	    if( p[a] )
		bracket_surface( p[a], px[a], px[(a+1)%4], xs[b++] );
	    else
		bracket_surface( p[(a+1)%4], px[(a+1)%4], px[a], xs[b++] );
	}
    }

    Vec3D n( xs[1][1]-xs[0][1], xs[0][0]-xs[1][0], 0.0 );
    if( !pol )
	n *= -1.0;

    n.normalize();

    //std::cout << "normal = " << n << "\n";

    return( n );
}


bool Geometry::vac_or_neu( int32_t i, int32_t j, int32_t k )
{
    signed char a = mesh_check( i, j, k );
    if( a >= -6 && a <= 0 ) // Vacuum or Neumann
	return( true );
    return( false );
}


void Geometry::build_mesh( void )
{
    int a;
    int32_t i, j, k;
    double x, y, z;

    if( ibsimu.get_verbose_output() ) {
	ibsimu.vout() << "Building mesh\n";
	ibsimu.vout().flush();
    }

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
	// Mark zmax
	if( _bound[5].type == BOUND_NEUMANN ) a = -6;
	else a = 6;
	for( j = 0; j < _size[1]; j++ ) {
	    for( i = 0; i < _size[0]; i++ ) {
		if( mesh(i,j,_size[2]-1) == 0 )
		    mesh(i,j,_size[2]-1) = a;
	    }
	}
    }

    // Mark solid edges
    for( k = 0; k < _size[2]; k++ ) {
	for( j = 0; j < _size[1]; j++ ) {
	    for( i = 0; i < _size[0]; i++ ) {
		// Check if solid node is an solid edge node
		if( (a = mesh(i,j,k)) >= 7 &&
		    ( vac_or_neu(i-1,j,  k  ) ||
		      vac_or_neu(i+1,j,  k  ) ||
		      vac_or_neu(i,  j-1,k  ) || 
		      vac_or_neu(i,  j+1,k  ) ||
		      vac_or_neu(i,  j,  k-1) || 
		      vac_or_neu(i,  j,  k+1) ) )  {
		    // Mark node as an edge
		    mesh(i,j,k) = -a;
		}
	    }
	}
    }

    // Report node counts
    if( ibsimu.get_verbose_output() ) {
	int b;
	int nc = nodecount();
	int vac = 0;
	int neu = 0;
	int dir = 0;
	int solid[_n];
	for( a = 0; a < (int)_n; a++ )
	    solid[a] = 0;

	for( a = 0; a < nc; a++ ) {
	    if( mesh(a) == 0 )
		vac++;
	    else if( (b = abs(mesh(a))) >= 7 )
		solid[b-7]++;
	    else if( mesh(a) < 0 )
		neu++;
	    else
		dir++;
	}
	
	ibsimu.vout() << "  Done. Built mesh with:\n";
	ibsimu.vout() << "  " << vac << " vacuum nodes\n";
	ibsimu.vout() << "  " << neu << " neumann nodes\n";
	ibsimu.vout() << "  " << dir << " dirichlet nodes\n";
	for( a = 0; a < (int)_n; a++ )
	    ibsimu.vout() << "  " << solid[a] << " solid " << a+7 << " nodes\n";
    }
}


void Geometry::save( const std::string &filename ) const
{
    if( ibsimu.get_verbose_output() )
	ibsimu.vout() << "Saving Geometry to file \'" << filename << "\'.\n";

    std::ofstream os( filename.c_str() );
    if( !os.good() )
	throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + filename + "\' for writing" ) );
    save( os );
    os.close();
}


void Geometry::save( std::ostream &os ) const
{
    Mesh::save( os );
    write_int32( os, _n );
    for( uint32_t a = 0; a < _n; a++ )
	_sdata[a]->save( os );
    for( uint32_t a = 0; a < _n+6; a++ )
	_bound[a].save( os );

    write_int8( os, _built );
    write_compressed_block( os, _size[0]*_size[1]*_size[2]*sizeof(signed char), _smesh );
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
	os << "bound[" << a << "] = " << _bound[a] << "\n";
    }
    os << "built = " << _built << "\n";
    os << "smesh = \n";
    for( int32_t j = 0; j < _size[1]; j++ ) {
	for( int32_t i = 0; i < _size[0]; i++ )
	    os << std::setw(2) << (int)mesh(i,j,0) << " ";
	os << "\n";
    }
}


