/*! \file solidgraph.cpp
 *  \brief Source code for solidgraph.cpp
 */

/* Copyright (c) 2005-2010 Taneli Kalvas. All rights reserved.
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

#include <vector>
#include <limits>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include "compmath.hpp"
#include "solidgraph.hpp"
#include "vec3d.hpp"
#include "lineclip.hpp"
#include "ibsimu.hpp"


int SolidGraph::get_mesh( const int i[3], int offsetx, int offsety ) const
{
    int j[3] = { i[0], i[1], i[2] };

    j[_vb[0]] += offsetx;
    j[_vb[1]] += offsety;

    return( _g.mesh_check( j[0], j[1], j[2] ) );
}


/* Return true if node at (i,j,k) is solid number N and is an edge
 * node. Here a node is counted as an edge node if some of the four
 * nearest neighbours differs from solid number N. Therefore also
 * simulation area boundaries will be counted as edges. 
 */
bool SolidGraph::is_edge( int N, const int i[3] ) const
{
    return( abs(get_mesh( i, 0, 0 )) == N &&
	    ( abs(get_mesh( i, -1,  0 )) != N ||
	      abs(get_mesh( i, +1,  0 )) != N ||
	      abs(get_mesh( i,  0, -1 )) != N ||
	      abs(get_mesh( i,  0, +1 )) != N ) );
}


void SolidGraph::build_solid( SolidPoints *solid, const int j[3], char *done, int lastN, int N )
{
    int next[3]  = {j[0], j[1], j[2]};
    int i[3]     = {j[0], j[1], j[2]};
    int first[3] = {j[0], j[1], j[2]};
    int a, dir, save;
    int loop_done = 0;
    
    //std::cout << solid << "\n";

    // Initialize direction
    if( lastN != N )
	dir = 9; // Going into solid -> direction 6 first
    else
	dir = 3; // Going out of solid -> direction 0 first

    //std::cout << "  Starting with from ("
    // << i[0] << ","
    // << i[1] << ","
    // << i[2] << ") with dir= " << dir << "\n";
    while( 1 ) {

	// Mark node done
	done[ i[_vb[0]] + i[_vb[1]]*_g.size(_vb[0]) ] = 1;

	/* One step counterclockwise from the direction we came from */
	dir -= 3;
	if( dir < 0 )
	    dir += 8;
	
	/* Go through directions */
	for( a = 0; a < 8; a++ ) {

	    //std::cout << "  Testing dir= " << dir << "\n";
	    save = 0;
	    switch( dir ) {
	    case 0:
		save = 1;
		next[_vb[0]] = i[_vb[0]];
		next[_vb[1]] = i[_vb[1]]+1;
		break;
	    case 1:
		if( abs(get_mesh( i,  0, +1 )) != N &&
		    abs(get_mesh( i, -1,  0 )) != N )
		    save = 1;
		next[_vb[0]] = i[_vb[0]]-1;
		next[_vb[1]] = i[_vb[1]]+1;
		break;
	    case 2:
		save = 1;
		next[_vb[0]] = i[_vb[0]]-1;
		next[_vb[1]] = i[_vb[1]];
		break;
	    case 3:
		if( abs(get_mesh( i,  0, -1 )) != N &&
		    abs(get_mesh( i, -1,  0 )) != N )
		    save = 1;
		next[_vb[0]] = i[_vb[0]]-1;
		next[_vb[1]] = i[_vb[1]]-1;
		break;
	    case 4:
		save = 1;
		next[_vb[0]] = i[_vb[0]];
		next[_vb[1]] = i[_vb[1]]-1;
		break;
	    case 5:
		if( abs(get_mesh( i,  0, -1 )) != N &&
		    abs(get_mesh( i, +1,  0 )) != N )
		    save = 1;
		next[_vb[0]] = i[_vb[0]]+1;
		next[_vb[1]] = i[_vb[1]]-1;
		break;
	    case 6:
		save = 1;
		next[_vb[0]] = i[_vb[0]]+1;
		next[_vb[1]] = i[_vb[1]];
		break;
	    case 7:
		if( abs(get_mesh( i,  0, +1 )) != N &&
		    abs(get_mesh( i, +1,  0 )) != N )
		    save = 1;
		next[_vb[0]] = i[_vb[0]]+1;
		next[_vb[1]] = i[_vb[1]]+1;
		break;
	    default:
		throw( Error( ERROR_LOCATION, "algorithm error" ) );
	    }

	    // If loop is done and direction same as starting direction
	    if( loop_done && ((dir == 6 && lastN != N) || (dir == 0 && lastN == N)) )
		break;

	    // If next point is an edge point, proceed to the next point without saving
	    if( abs(_g.mesh_check( next[0], next[1], next[2] )) == N )
		break;

	    // Save solid edge point between i and next
	    if( save ) {
		Vec3D x1( i[0]*_g.h()+_g.origo(0),
			  i[1]*_g.h()+_g.origo(1),
			  i[2]*_g.h()+_g.origo(2) );
		Vec3D x2( next[0]*_g.h()+_g.origo(0),
			  next[1]*_g.h()+_g.origo(1),
			  next[2]*_g.h()+_g.origo(2) );
		Vec3D xsurf;
		_g.bracket_surface( N, x1, x2, xsurf );
		//std::cout << "  Saving point (" 
		//	  << xsurf(vb[0]) << ","
		//	  << xsurf(vb[1]) << ")\n";
		//std::cout << solid->N << "\n";
		//solid->p.reserve( 10 );
		//std::cout << solid->p.size() << " " << solid->p.capacity() << "\n";
		solid->p.push_back( Point( xsurf(_vb[0]), xsurf(_vb[1]) ) );
		//std::cout << "  done\n";
	    }

	    // Progress direction counterclockwise
	    dir++;
	    if( dir > 7 )
		dir -= 8;
	}
	if( a == 8 || loop_done ) {
	    //std::cout << "  Loop done\n";
	    break;
	}

	// Advance to next point
	//std::cout << "  Advance to (" 
	//	  << i[0] << ","
	//	  << i[1] << ","
	//	  << i[2] << ")\n";
	i[0] = next[0];
	i[1] = next[1];
	i[2] = next[2];

	if( i[0] == first[0] && i[1] == first[1] && i[2] == first[2] ) {
	    loop_done = 1;
	}
    }
}


void SolidGraph::clear_data( void )
{
    for( size_t a = 0; a < _solid.size(); a++  ) {
	delete _solid[a];
    }
    _solid.clear();
}


void SolidGraph::build_data( void )
{
    // Clear old plot
    clear_data();

    char *done = new char[ _g.size(_vb[0]) * _g.size(_vb[1]) ];
    memset( done, 0, _g.size(_vb[0])*_g.size(_vb[1])*sizeof(char) );

    // Go through mesh
    int lastN = 0;
    int i[3];
    i[_vb[2]] = (int)floor(_level+0.5);
    for( i[_vb[1]] = 0; i[_vb[1]] < _g.size(_vb[1]); i[_vb[1]]++ ) {
	for( i[_vb[0]] = 0; i[_vb[0]] < _g.size(_vb[0]); i[_vb[0]]++ ) {

	    //std::cout << "Processing (" 
	    // << i[0] << "," 
	    // << i[1] << "," 
	    // << i[2] << "): ";

	    int N = abs( _g.mesh( i[0], i[1], i[2] ) );
	    //std::cout << "N = " << N << " ";
	    // Skip processed nodes and vacuum
	    if( N < 7 ) {
		//std::cout << "vacuum or boundary\n";
		continue;
	    } else  if( done[ i[_vb[0]] + i[_vb[1]]*_g.size(_vb[0]) ] ) {
		//std::cout << "done\n";
		continue;
	    } else if( !is_edge( N, i ) ) {
		//std::cout << "not an edge\n";
		continue;
	    }

	    // If solid is a continuation of processed solid
	    size_t a;
	    for( a = 0; a < _solid.size(); a++ ) {
		if( N == _solid[a]->N ) {
		    // Add break
		    //std::cout << "building continuation\n";
		    _solid[a]->p.push_back( Point(std::numeric_limits<double>::quiet_NaN(),
						  std::numeric_limits<double>::quiet_NaN()) );
		    build_solid( _solid[a], i, done, lastN, N );
		    break;
		}
	    }
	    if( a == _solid.size() ) {
		// Create a new solid
		//std::cout << "building new\n";
		_solid.push_back( new SolidPoints( N ) );
		//std::cout << "push_back done\n";
		//std::cout << solid[a] << "\n";
		build_solid( _solid[a], i, done, lastN, N );
	    }
	    lastN = N;

	}
    }

    delete [] done;
}




SolidGraph::SolidGraph( const Geometry &g ) 
  : _color(Color(0.2,0.2,1.0)), _g(g), _cache(true)
{
}


SolidGraph::~SolidGraph()
{
    for( size_t a = 0; a < _solid.size(); a++ )
	delete _solid[a];
}


void SolidGraph::disable_cache( void )
{
    _cache = false;
}


void SolidGraph::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] )
{
    //if( ibsimu.get_verbose_output() )
    //std::cout << "  Plotting solids\n";

    //std::cout << "--Plotting solids--\n";
    if( !_cache || _solid.size() == 0 || _oview != _view || _olevel != _level ) {
	// First round or change happened
	//_ilevel = (int)floor( (_level - _g.origo(_vb[2])) / _g.h() + 0.5 );
	//std::cout << "--Building solid plot data for level " << _level << "--\n";
	build_data();
    }
    _oview = _view;
    _olevel = _level;

    // Set drawing properties
    cairo_set_source_rgba( cairo, _color[0], _color[1], _color[2], _color[3] );
    cairo_set_line_width( cairo, 1.0 );

    // Set clipping ranges
    double clip[4];
    cm->transform( &clip[0], &range[0] );
    cm->transform( &clip[2], &range[2] );
    LineClip lc( cairo );
    lc.set( clip[0], clip[1], clip[2], clip[3] );

    for( size_t a = 0; a < _solid.size(); a++ ) {

	SolidPoints *sp = _solid[a];
	double xout[2];
	cm->transform( xout, sp->p[0].x );
	lc.move_to( xout[0], xout[1] );

	for( size_t b = 1; b < sp->p.size(); b++ ) {
	    if( comp_isnan( sp->p[b].x[0] ) ) {
		// Break in path, separate paths
		do b++;
		while( b != sp->p.size() && comp_isnan( sp->p[b].x[0] ) );
		if( b == sp->p.size() )
		    break;
		lc.close_path();
		cm->transform( xout, sp->p[b].x );
		lc.move_to( xout[0], xout[1] );
	    } else {
		cm->transform( xout, sp->p[b].x );
		lc.line_to( xout[0], xout[1] );
	    }
	}

	lc.fill();
    }

}


void SolidGraph::get_bbox( double bbox[4] )
{
    bbox[0] = _g.origo( _vb[0] );
    bbox[1] = _g.origo( _vb[1] );
    bbox[2] = _g.max( _vb[0] );
    bbox[3] = _g.max( _vb[1] );
}



















