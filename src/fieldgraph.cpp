/*! \file fieldgraph.cpp
 *  \brief Source code for fieldgraph.cpp
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

#include "fieldgraph.hpp"
#include "ibsimu.hpp"


FieldGraph::FieldGraph( const ScalarField *field )
    : _scalarfield(field), _colormap(NULL), _enabled(true)
{
    
}


FieldGraph::~FieldGraph()
{
    if( _colormap )
	delete _colormap;
}


void FieldGraph::enable( bool enable )
{
    _enabled = enable;
    if( !_enabled && _colormap )
	delete _colormap;
}


void FieldGraph::build_scalarfield_plot( void )
{
    if( _scalarfield ) {
	double range[4] = { _scalarfield->origo( _vb[0] ),
			    _scalarfield->origo( _vb[1] ),
			    _scalarfield->max  ( _vb[0] ),
			    _scalarfield->max  ( _vb[1] ) };
	size_t n = _scalarfield->size( _vb[0] );
	size_t m = _scalarfield->size( _vb[1] );
	std::vector<double> data;
	data.reserve( n*m );
	size_t ind[3];
	ind[_vb[2]] = (size_t)floor(_level+0.5);
	for( size_t j = 0; j < m; j++ ) {
	    ind[_vb[1]] = j;
	    for( size_t i = 0; i < n; i++ ) {
		ind[_vb[0]] = i;
		data.push_back( (*_scalarfield)( ind[0], ind[1], ind[2] ) );
	    }
	}
	_colormap = new Colormap( range, n, m, data );
	std::vector<Palette::Entry> pentry;
	pentry.push_back( Palette::Entry( Color(1,1,1), 0 ) );
	pentry.push_back( Palette::Entry( Color(1,1,0), 1 ) );
	pentry.push_back( Palette::Entry( Color(1,0,0), 2 ) );
	pentry.push_back( Palette::Entry( Color(0,0,0), 3 ) );
	Palette p( pentry );
	_colormap->set_palette( p );
	_colormap->set_zscale( ZSCALE_RELLOG );
    }
}


void FieldGraph::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] )
{
    if( ibsimu.get_verbose_output() )
	std::cout << "  Plotting field\n";

    if( _colormap == NULL || _oview != _view || _olevel != _level ) {
	// First plot or changed happened
	if( _scalarfield && _enabled )
	    build_scalarfield_plot();
    }
    _oview = _view;
    _olevel = _level;

    if( _colormap )
	_colormap->plot( cairo, cm, range );
}


void FieldGraph::get_bbox( double bbox[4] )
{
    if( _scalarfield ) {
	bbox[0] = _scalarfield->origo( _vb[0] );
	bbox[1] = _scalarfield->origo( _vb[1] );
	bbox[2] = _scalarfield->max( _vb[0] );
	bbox[3] = _scalarfield->max( _vb[1] );
    }
}


















