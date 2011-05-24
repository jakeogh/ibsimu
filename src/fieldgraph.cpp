/*! \file fieldgraph.cpp
 *  \brief %Graph for plotting fields
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

#include "fieldgraph.hpp"
#include "ibsimu.hpp"


FieldGraph::FieldGraph( const ScalarField *field )
    : _field_type(FIELD_NONE), _geom(NULL), _scalarfield(field), _vectorfield(NULL), 
      _colormap(NULL), _enabled(true)
{
    
}


FieldGraph::FieldGraph( const Geometry *geom, const VectorField *field, field_type_e field_type )
    : _field_type(field_type), _geom(geom), _scalarfield(NULL), _vectorfield(field), 
      _colormap(NULL), _enabled(true)
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


void FieldGraph::build_vectorfield_plot( void )
{
    if( !_vectorfield )
	return;

    // Build data for colormap based on geometry size
    double range[4] = { _geom->origo( _vb[0] ),
			_geom->origo( _vb[1] ),
			_geom->max  ( _vb[0] ),
			_geom->max  ( _vb[1] ) };
    size_t n = _geom->size( _vb[0] );
    size_t m = _geom->size( _vb[1] );

    std::vector<double> data;
    data.reserve( n*m );
    //size_t ind[3];
    Vec3D x;
    //ind[_vb[2]] = _level;
    x[_vb[2]] = _geom->origo(_vb[2]) + _level*_geom->h();
    for( size_t j = 0; j < m; j++ ) {
	//ind[_vb[1]] = j;
	x[_vb[1]] = _geom->origo(_vb[1]) + j*_geom->h();
	for( size_t i = 0; i < n; i++ ) {
	    //ind[_vb[0]] = i;
	    x[_vb[0]] = _geom->origo(_vb[0]) + i*_geom->h();

	    // Make data for point (ind[0], ind[1], ind[2])
	    //data.push_back( (*_scalarfield)( ind[0], ind[1], ind[2] ) );

	    Vec3D F = (*_vectorfield)( x );

	    switch( _field_type ) {
	    case FIELD_EFIELD:
	    case FIELD_BFIELD:
		data.push_back( F.norm2() );
		break;
	    case FIELD_EFIELD_X:
	    case FIELD_BFIELD_X:
		data.push_back( F[0] );
		break;
	    case FIELD_EFIELD_Y:
	    case FIELD_BFIELD_Y:
		data.push_back( F[1] );
		break;
	    case FIELD_EFIELD_Z:
	    case FIELD_BFIELD_Z:
		data.push_back( F[2] );
		break;
	    default:
		throw( Error( ERROR_LOCATION, "unknown field type" ) );
		break;
	    }
	}
    }

    // Set up colormap
    _colormap = new Colormap( range, n, m, data );
    //double zmin, zmax;
    //_scalarfield->get_minmax( zmin, zmax );
    //_colormap->set_zrange( zmin, zmax );
    //double zspan = zmax - zmin;

    // Set palette
    std::vector<Palette::Entry> pentry;
    pentry.push_back( Palette::Entry( Color(0,0,0), -1.00 ) );
    pentry.push_back( Palette::Entry( Color(0,0,1), -0.67 ) );
    pentry.push_back( Palette::Entry( Color(0,1,1), -0.33 ) );
    pentry.push_back( Palette::Entry( Color(1,1,1),  0.00 ) );
    pentry.push_back( Palette::Entry( Color(1,1,0),  0.33 ) );
    pentry.push_back( Palette::Entry( Color(1,0,0),  0.67 ) );
    pentry.push_back( Palette::Entry( Color(0,0,0),  1.00 ) );
    Palette p( pentry );
    _colormap->set_palette( p );
    if( _logscale )
	_colormap->set_zscale( ZSCALE_RELLOG );
}


void FieldGraph::build_scalarfield_plot( void )
{
    if( !_scalarfield )
	return;

    // Build data for colormap based on scalarfield size
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

    // Set up colormap
    _colormap = new Colormap( range, n, m, data );
    double zmin, zmax;
    _scalarfield->get_minmax( zmin, zmax );
    _colormap->set_zrange( zmin, zmax );
    double zspan = zmax - zmin;
    
    //std::cout << "Colormap zmin = " << zmin << "\n";
    //std::cout << "Colormap zmax = " << zmax << "\n";
    //std::cout << "Colormap zspan = " << zspan << "\n";

    // Set palette
    std::vector<Palette::Entry> pentry;
    if( zmin >= -1.0e-6*zspan && zmax >= 0.0 ) {
	// Palette for positive beam
	//std::cout << "Palette for positive beam\n";
	pentry.push_back( Palette::Entry( Color(1,1,1), 0 ) );
	pentry.push_back( Palette::Entry( Color(1,1,0), 1 ) );
	pentry.push_back( Palette::Entry( Color(1,0,0), 2 ) );
	pentry.push_back( Palette::Entry( Color(0,0,0), 3 ) );
    } else if( zmax <= 1.0e-6*zspan && zmin <= 0.0 ) {
	// Palette for negative beam
	//std::cout << "Palette for negative beam\n";
	pentry.push_back( Palette::Entry( Color(1,1,1), 3 ) );
	pentry.push_back( Palette::Entry( Color(1,1,0), 2 ) );
	pentry.push_back( Palette::Entry( Color(1,0,0), 1 ) );
	pentry.push_back( Palette::Entry( Color(0,0,0), 0 ) );
    } else {
	// Palette for positive and negative beam
	//std::cout << "Palette for dual polarity beam\n";
	pentry.push_back( Palette::Entry( Color(0,0,0), zmin ) );
	pentry.push_back( Palette::Entry( Color(0,0,1), 0.67*zmin ) );
	pentry.push_back( Palette::Entry( Color(0,1,1), 0.33*zmin ) );
	pentry.push_back( Palette::Entry( Color(1,1,1), 0 ) );
	pentry.push_back( Palette::Entry( Color(1,1,0), 0.33*zmax ) );
	pentry.push_back( Palette::Entry( Color(1,0,0), 0.67*zmax ) );
	pentry.push_back( Palette::Entry( Color(0,0,0), zmax ) );
    }
    Palette p( pentry );
    _colormap->set_palette( p );
    if( _logscale )
	_colormap->set_zscale( ZSCALE_RELLOG );
}


void FieldGraph::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] )
{
    if( _colormap == NULL || _oview != _view || _olevel != _level ) {
	// First plot or changed view happened

	// Scalarfield
	if( _scalarfield && _enabled )
	    build_scalarfield_plot();

	// Vectorfield
	else if( _vectorfield && _enabled )
	    build_vectorfield_plot();
    }
    _oview = _view;
    _olevel = _level;

    if( _colormap )
	_colormap->plot( cairo, cm, range );
}


void FieldGraph::plot_sample( cairo_t *cairo, double x[2] )
{

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
