/*! \file solidgraph.cpp
 *  \brief %Graph for plotting solids.
 */

/* Copyright (c) 2005-2012 Taneli Kalvas. All rights reserved.
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


//#define DEBUG_SOLIDGRAPH 1


uint32_t SolidGraph::get_mesh( const int32_t i[3], int offsetx, int offsety ) const
{
    int j[3] = { i[0], i[1], i[2] };

    j[_vb[0]] += offsetx;
    j[_vb[1]] += offsety;

    return( _geom.mesh_check( j[0], j[1], j[2] ) );
}


/* Return true if node at (i,j,k) is solid mesh node number N and is
 * an edge node. Here a node is counted as an edge node if some of the
 * four nearest neighbours differs from solid number N. Therefore also
 * simulation area boundaries will be counted as edges.
 */
bool SolidGraph::is_edge( uint32_t node, const int32_t i[3] ) const
{
    return( get_mesh( i, 0, 0 ) == node &&
	    ( get_mesh( i, -1,  0 ) != node ||
	      get_mesh( i, +1,  0 ) != node ||
	      get_mesh( i,  0, -1 ) != node ||
	      get_mesh( i,  0, +1 ) != node ) );
}


/* Loop around solid and save points to \a solid. Uses near solid data
 * from Geometry for solid surface location when available, otherwise
 * draws on the solid edge node. Only does drawing in coordinate axes
 * directions from the surface. Direction of travel includes 45 degree
 * directions:
 *
 *   1 0 7
 *   2 X 6
 *   3 4 5
 */
void SolidGraph::build_solid( SolidPoints *solid, const int32_t j[3], char *done, bool out, uint32_t node )
{
    int32_t next[3]  = {j[0], j[1], j[2]};
    int32_t i[3]     = {j[0], j[1], j[2]};
    int32_t first[3] = {j[0], j[1], j[2]};
    int a, dir;
    int loop_done = 0;

#ifdef DEBUG_SOLIDGRAPH    
    std::cout << "\nbuild_solid()\n";
#endif

    // Initialize direction
    if( out )
	dir = 3; // Going out of solid -> direction 0 first (X will be subtracted later)
    else
	dir = 9; // Going into solid -> direction 3 first (X will be subtracted later)

    while( 1 ) {

	// Mark node done
	done[ i[_vb[0]] + i[_vb[1]]*_geom.size(_vb[0]) ] = 1;

	// One step counterclockwise from the direction we came from
	dir -= 3;
	if( dir < 0 )
	    dir += 8;
	
	// Go through directions
	for( a = 0; a < 8; a++ ) {

#ifdef DEBUG_SOLIDGRAPH
	    std::cout << "  Testing dir = " << dir << "\n";
#endif

	    switch( dir ) {
	    case 0:
		next[_vb[0]] = i[_vb[0]];
		next[_vb[1]] = i[_vb[1]]+1;
		break;
	    case 1:
		next[_vb[0]] = i[_vb[0]]-1;
		next[_vb[1]] = i[_vb[1]]+1;
		break;
	    case 2:
		next[_vb[0]] = i[_vb[0]]-1;
		next[_vb[1]] = i[_vb[1]];
		break;
	    case 3:
		next[_vb[0]] = i[_vb[0]]-1;
		next[_vb[1]] = i[_vb[1]]-1;
		break;
	    case 4:
		next[_vb[0]] = i[_vb[0]];
		next[_vb[1]] = i[_vb[1]]-1;
		break;
	    case 5:
		next[_vb[0]] = i[_vb[0]]+1;
		next[_vb[1]] = i[_vb[1]]-1;
		break;
	    case 6:
		next[_vb[0]] = i[_vb[0]]+1;
		next[_vb[1]] = i[_vb[1]];
		break;
	    case 7:
		next[_vb[0]] = i[_vb[0]]+1;
		next[_vb[1]] = i[_vb[1]]+1;
		break;
	    }

#ifdef DEBUG_SOLIDGRAPH
	    std::cout << "    save = " << save << "\n";
	    std::cout << "    next = (" 
		      << next[0] << ", "
		      << next[1] << ", "
		      << next[2] << ")\n";
#endif

	    // If loop is done and direction same as starting direction
	    if( loop_done && ((dir == 6 && !out) || (dir == 0 && out)) )
		break;

	    // If next point is an edge point, proceed to the next
	    // point without saving a point
	    if( is_edge( node, next ) )
		break;

	    // Save solid surface location between i and next
	    if( dir == 0 || dir == 2 || dir == 4 || dir == 6 ) {

		uint32_t mesh = _geom.mesh_check( next[0], next[1], next[2] );
		if( (mesh & SMESH_NODE_ID_MASK) == SMESH_NODE_ID_NEAR_SOLID ) {
		    uint32_t dir_to_solid;
		    if( dir == 0 ) {
			dir_to_solid = 2*_vb[1];
			double dist = (255-_geom.solid_dist( next[0], next[1], next[2], dir_to_solid ))/255.0;
			solid->p.push_back( Point( i[_vb[0]]*_geom.h()+_geom.origo(_vb[0]),
						   (i[_vb[1]]+dist)*_geom.h()+_geom.origo(_vb[1]) ) );
		    } else if( dir == 2 ) {
			dir_to_solid = 2*_vb[0]+1;
			double dist = (255-_geom.solid_dist( next[0], next[1], next[2], dir_to_solid ))/255.0;
			solid->p.push_back( Point( (i[_vb[0]]-dist)*_geom.h()+_geom.origo(_vb[0]),
						   i[_vb[1]]*_geom.h()+_geom.origo(_vb[1]) ) );
		    } else if( dir == 4 ) {
			dir_to_solid = 2*_vb[1]+1;
			double dist = (255-_geom.solid_dist( next[0], next[1], next[2], dir_to_solid ))/255.0;
			solid->p.push_back( Point( i[_vb[0]]*_geom.h()+_geom.origo(_vb[0]),
						   (i[_vb[1]]-dist)*_geom.h()+_geom.origo(_vb[1]) ) );
		    } else { // dir == 6
			dir_to_solid = 2*_vb[0];
			double dist = (255-_geom.solid_dist( next[0], next[1], next[2], dir_to_solid ))/255.0;
			solid->p.push_back( Point( (i[_vb[0]]+dist)*_geom.h()+_geom.origo(_vb[0]),
						   i[_vb[1]]*_geom.h()+_geom.origo(_vb[1]) ) );
		    }
		} else {
		    // Save edge node location (would halfway work better?)
		    solid->p.push_back( Point( i[_vb[0]]*_geom.h()+_geom.origo(_vb[0]),
					       i[_vb[1]]*_geom.h()+_geom.origo(_vb[1]) ) );
		}
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

    // Array for marking processed nodes
    char *done = new char[ _geom.size(_vb[0]) * _geom.size(_vb[1]) ];
    memset( done, 0, _geom.size(_vb[0])*_geom.size(_vb[1])*sizeof(char) );

    // Go through mesh
    int32_t i[3];
    uint32_t lastnode;
    uint32_t node = SMESH_NODE_ID_PURE_VACUUM;
    i[_vb[2]] = (int)floor(_level+0.5);
    for( i[_vb[1]] = 0; i[_vb[1]] < (int32_t)_geom.size(_vb[1]); i[_vb[1]]++ ) {
	for( i[_vb[0]] = 0; i[_vb[0]] < (int32_t)_geom.size(_vb[0]); i[_vb[0]]++ ) {

	    //std::cout << "Processing (" 
	    // << i[0] << "," 
	    // << i[1] << "," 
	    // << i[2] << "): ";

	    // Skip until unprocessed edge is found
	    lastnode = node;
	    node = _geom.mesh( i[0], i[1], i[2] );
	    uint32_t nodeid = node & SMESH_NODE_ID_MASK;
	    if( nodeid == SMESH_NODE_ID_NEAR_SOLID ||
		nodeid == SMESH_NODE_ID_NEAR_SOLID_FIX ||
		nodeid == SMESH_NODE_ID_PURE_VACUUM ||
		nodeid == SMESH_NODE_ID_PURE_VACUUM_FIX ||
		nodeid == SMESH_NODE_ID_NEUMANN ) {
		//std::cout << "vacuum\n";
		continue;
	    } else if( nodeid == SMESH_NODE_ID_DIRICHLET &&
		       (node & SMESH_BOUNDARY_NUMBER_MASK) < 7 ) {
		//std::cout << "Dirichlet simulation box boundary\n";
		continue;
	    } else if( done[ i[_vb[0]] + i[_vb[1]]*_geom.size(_vb[0]) ] ) {
		//std::cout << "done\n";
		continue;
	    } else if( !is_edge( node, i ) ) {
		//std::cout << "not an edge\n";
		continue;
	    }

	    // If solid is a continuation of processed solid
	    size_t a;
	    for( a = 0; a < _solid.size(); a++ ) {
		if( (node & SMESH_BOUNDARY_NUMBER_MASK) == _solid[a]->N ) {
		    // Add break
		    //std::cout << "building continuation\n";
		    _solid[a]->p.push_back( Point(std::numeric_limits<double>::quiet_NaN(),
						  std::numeric_limits<double>::quiet_NaN()) );
		    build_solid( _solid[a], i, done, (lastnode == node), node );
		    break;
		}
	    }
	    if( a == _solid.size() ) {
		// Create a new solid
		//std::cout << "building new\n";
		_solid.push_back( new SolidPoints( node & SMESH_BOUNDARY_NUMBER_MASK ) );
		//std::cout << "push_back done\n";
		//std::cout << solid[a] << "\n";
		//build_solid( _solid[a], i, done, (lastnode == node), node );
		build_solid( _solid[a], i, done, (lastnode == node), node );
	    }
	}
    }

    delete [] done;
}


SolidGraph::SolidGraph( const Geometry &geom ) 
    : Graph3D(geom), _geom(geom), _color(Color(0.2,0.2,1.0)), _cache(true)
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
    if( !_cache || _solid.size() == 0 || _oview != _view || _olevel != _level ) {
	// First round or change happened
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


void SolidGraph::plot_sample( cairo_t *cairo, double x, double y, double width, double height )
{

}


void SolidGraph::get_bbox( double bbox[4] )
{
    bbox[0] = _geom.origo( _vb[0] );
    bbox[1] = _geom.origo( _vb[1] );
    bbox[2] = _geom.max( _vb[0] );
    bbox[3] = _geom.max( _vb[1] );
}
