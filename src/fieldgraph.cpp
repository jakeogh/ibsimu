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

#include <limits>
#include "fieldgraph.hpp"
#include "meshvectorfield.hpp"
#include "ibsimu.hpp"


FieldGraph::FieldGraph( const Geometry &geom )
    : Graph3D(geom), _geom(geom), _field_type(FIELD_NONE), _scalarfield(NULL),
      _vectorfield(NULL), _first(true), _enabled(true), _zmin(0.0), _zmax(0.0)
{
    
}


FieldGraph::FieldGraph( const Geometry &geom, field_type_e field_type, const ScalarField *field )
    : Graph3D(geom), _geom(geom), _first(true), _enabled(true)
{
    set_field( field_type, field );
}


FieldGraph::FieldGraph( const Geometry &geom, field_type_e field_type, const VectorField *field )
    : Graph3D(geom), _geom(geom), _first(true), _enabled(true)
{
    set_field( field_type, field );
}


FieldGraph::~FieldGraph()
{
}


field_type_e FieldGraph::field_type( void )
{
    return( _field_type );
}


void FieldGraph::set_field( field_type_e field_type, const ScalarField *field )
{
    // Set parameters
    _first = true;
    _field_type = field_type;
    _scalarfield = field;
    _vectorfield = NULL;

    if( field == NULL ) {
	_field_type = FIELD_NONE;
	return;
    }

    _scalarfield->get_minmax( _zmin, _zmax );

    // Set default colormap palette
    std::vector<Palette::Entry> pentry;
    double zspan = _zmax - _zmin;
    if( _zmin >= -1.0e-6*zspan && _zmax >= 0.0 ) {
	// Red palette
	pentry.push_back( Palette::Entry( Color(1,1,1), 0 ) );
	pentry.push_back( Palette::Entry( Color(1,1,0), 1 ) );
	pentry.push_back( Palette::Entry( Color(1,0,0), 2 ) );
	pentry.push_back( Palette::Entry( Color(0,0,0), 3 ) );
    } else if( _zmax <= 1.0e-6*zspan && _zmin <= 0.0 ) {
	// Red palette
	pentry.push_back( Palette::Entry( Color(1,1,1), 3 ) );
	pentry.push_back( Palette::Entry( Color(1,1,0), 2 ) );
	pentry.push_back( Palette::Entry( Color(1,0,0), 1 ) );
	pentry.push_back( Palette::Entry( Color(0,0,0), 0 ) );
    } else {
	// Mixed palette, forcing zero at white
	pentry.push_back( Palette::Entry( Color(0,0,0), _zmin ) );
	pentry.push_back( Palette::Entry( Color(0,0,1), 0.67*_zmin ) );
	pentry.push_back( Palette::Entry( Color(0,1,1), 0.33*_zmin ) );
	pentry.push_back( Palette::Entry( Color(1,1,1), 0 ) );
	pentry.push_back( Palette::Entry( Color(1,1,0), 0.33*_zmax ) );
	pentry.push_back( Palette::Entry( Color(1,0,0), 0.67*_zmax ) );
	pentry.push_back( Palette::Entry( Color(0,0,0), _zmax ) );
    }
    Palette p( pentry );
    set_palette( p );
    Colormap::set_zrange( _zmin, _zmax );
}


void FieldGraph::set_field( field_type_e field_type, const VectorField *field )
{
    // Set parameters
    _first = true;
    _field_type = field_type;
    _scalarfield = NULL;
    _vectorfield = field;

    if( field == NULL ) {
	_field_type = FIELD_NONE;
	return;
    }

    const MeshVectorField *mvf = dynamic_cast<const MeshVectorField *>( field );
    if( mvf ) {

	Vec3D min, max;
	switch( _field_type ) {
	case FIELD_EFIELD:
	case FIELD_BFIELD:
	    mvf->get_minmax( _zmin, _zmax );
	    break;
	case FIELD_EFIELD_X:
	case FIELD_BFIELD_X:
	    mvf->get_minmax( min, max );
	    _zmin = min[0];
	    _zmax = max[0];
	    break;
	case FIELD_EFIELD_Y:
	case FIELD_BFIELD_Y:
	    mvf->get_minmax( min, max );
	    _zmin = min[1];
	    _zmax = max[1];
	    break;
	case FIELD_EFIELD_Z:
	case FIELD_BFIELD_Z:
	    mvf->get_minmax( min, max );
	    _zmin = min[2];
	    _zmax = max[2];
	    break;
	default:
	    throw( Error( ERROR_LOCATION, "unknown field type" ) );
	    break;
	}

    } else {

	Vec3D min, max;
	switch( _field_type ) {
	case FIELD_EFIELD:
	case FIELD_BFIELD:
	    field->get_minmax( _geom, _zmin, _zmax );
	    break;
	case FIELD_EFIELD_X:
	case FIELD_BFIELD_X:
	    field->get_minmax( _geom, min, max );
	    _zmin = min[0];
	    _zmax = max[0];
	    break;
	case FIELD_EFIELD_Y:
	case FIELD_BFIELD_Y:
	    field->get_minmax( _geom, min, max );
	    _zmin = min[1];
	    _zmax = max[1];
	    break;
	case FIELD_EFIELD_Z:
	case FIELD_BFIELD_Z:
	    field->get_minmax( _geom, min, max );
	    _zmin = min[2];
	    _zmax = max[2];
	    break;
	default:
	    throw( Error( ERROR_LOCATION, "unknown field type" ) );
	    break;
	}

    }

    // Set default colormap palette
    std::vector<Palette::Entry> pentry;
    pentry.push_back( Palette::Entry( Color(0,0,0), -1.00 ) );
    pentry.push_back( Palette::Entry( Color(0,0,1), -0.67 ) );
    pentry.push_back( Palette::Entry( Color(0,1,1), -0.33 ) );
    pentry.push_back( Palette::Entry( Color(1,1,1),  0.00 ) );
    pentry.push_back( Palette::Entry( Color(1,1,0),  0.33 ) );
    pentry.push_back( Palette::Entry( Color(1,0,0),  0.67 ) );
    pentry.push_back( Palette::Entry( Color(0,0,0),  1.00 ) );
    Palette p( pentry );
    set_palette( p );
    Colormap::set_zrange( _zmin, _zmax );
}


void FieldGraph::enable( bool enable )
{
    _enabled = enable;
}


void FieldGraph::set_zrange( double min, double max )
{
    _zmin = min;
    _zmax = max;

    Colormap::set_zrange( _zmin, _zmax );
}


void FieldGraph::build_vectorfield_plot( void )
{
    // Build data for colormap based on geometry size
    double range[4] = { _geom.origo( _vb[0] ),
			_geom.origo( _vb[1] ),
			_geom.max  ( _vb[0] ),
			_geom.max  ( _vb[1] ) };
    size_t n = _geom.size( _vb[0] );
    size_t m = _geom.size( _vb[1] );

    std::vector<double> data;
    data.reserve( n*m );
    Vec3D x;
    x[_vb[2]] = _geom.origo(_vb[2]) + _level*_geom.h();
    for( size_t j = 0; j < m; j++ ) {
	x[_vb[1]] = _geom.origo(_vb[1]) + j*_geom.h();
	for( size_t i = 0; i < n; i++ ) {
	    x[_vb[0]] = _geom.origo(_vb[0]) + i*_geom.h();

	    Vec3D F = (*_vectorfield)( x );
	    double Fv = 0.0;

	    switch( _field_type ) {
	    case FIELD_EFIELD:
	    case FIELD_BFIELD:
		Fv = F.norm2();
		break;
	    case FIELD_EFIELD_X:
	    case FIELD_BFIELD_X:
		Fv = F[0];
		break;
	    case FIELD_EFIELD_Y:
	    case FIELD_BFIELD_Y:
		Fv = F[1];
		break;
	    case FIELD_EFIELD_Z:
	    case FIELD_BFIELD_Z:
		Fv = F[2];
		break;
	    default:
		throw( Error( ERROR_LOCATION, "unknown field type" ) );
		break;
	    }

	    data.push_back( Fv );
	}
    }

    // Set up colormap
    set_data( range, n, m, data );
    Colormap::set_zrange( _zmin, _zmax );
}


void FieldGraph::build_meshvectorfield_plot( void )
{
    const MeshVectorField *mvf = dynamic_cast<const MeshVectorField *>( _vectorfield );

    if( !mvf )
	return;

    // Build data for colormap based on scalarfield size
    double range[4] = { mvf->origo( _vb[0] ),
			mvf->origo( _vb[1] ),
			mvf->max  ( _vb[0] ),
			mvf->max  ( _vb[1] ) };
    size_t n = mvf->size( _vb[0] );
    size_t m = mvf->size( _vb[1] );

    std::vector<double> data;
    data.reserve( n*m );
    int32_t ind[3];
    ind[_vb[2]] = (int32_t)floor( (_level_si - mvf->origo(_vb[2])) / mvf->h() + 0.5 );

    // Check that field is visible
    if( ind[_vb[2]] < 0 || ind[_vb[2]] >= mvf->size( _vb[2] ) ) {
	clear_data();
	return;
    }

    for( size_t j = 0; j < m; j++ ) {
	ind[_vb[1]] = j;
	for( size_t i = 0; i < n; i++ ) {
	    ind[_vb[0]] = i;
	    Vec3D F = (*mvf)( ind[0], ind[1], ind[2] );
	    double Fv = 0.0;

	    switch( _field_type ) {
	    case FIELD_EFIELD:
	    case FIELD_BFIELD:
		Fv = F.norm2();
		break;
	    case FIELD_EFIELD_X:
	    case FIELD_BFIELD_X:
		Fv = F[0];
		break;
	    case FIELD_EFIELD_Y:
	    case FIELD_BFIELD_Y:
		Fv = F[1];
		break;
	    case FIELD_EFIELD_Z:
	    case FIELD_BFIELD_Z:
		Fv = F[2];
		break;
	    default:
		throw( Error( ERROR_LOCATION, "unknown field type" ) );
		break;
	    }

	    data.push_back( Fv );
	}
    }

    // Set up colormap
    set_data( range, n, m, data );
    Colormap::set_zrange( _zmin, _zmax );
}


void FieldGraph::build_scalarfield_plot( void )
{
    // Build data for colormap based on scalarfield size
    double range[4] = { _scalarfield->origo( _vb[0] ),
			_scalarfield->origo( _vb[1] ),
			_scalarfield->max  ( _vb[0] ),
			_scalarfield->max  ( _vb[1] ) };
    size_t n = _scalarfield->size( _vb[0] );
    size_t m = _scalarfield->size( _vb[1] );

    std::vector<double> data;
    data.reserve( n*m );
    int32_t ind[3];
    ind[_vb[2]] = (int32_t)floor( (_level_si - _scalarfield->origo(_vb[2])) / 
				  _scalarfield->h() + 0.5 );
    
    // Check that field is visible
    if( ind[_vb[2]] < 0 || ind[_vb[2]] >= _scalarfield->size( _vb[2] ) )
	clear_data();

    for( size_t j = 0; j < m; j++ ) {
	ind[_vb[1]] = j;
	for( size_t i = 0; i < n; i++ ) {
	    ind[_vb[0]] = i;
	    double Fv = (*_scalarfield)( ind[0], ind[1], ind[2] );
	    data.push_back( Fv );
	}
    }

    // Set up colormap
    set_data( range, n, m, data );
    Colormap::set_zrange( _zmin, _zmax );
}


void FieldGraph::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] )
{
    if( _first || _oview != _view || _olevel != _level ) {
	// First plot or changed view happened

	// Scalarfield
	if( _scalarfield && _enabled )
	    build_scalarfield_plot();

	else if( _vectorfield && _enabled ) {
	    //const MeshVectorField *mvf = dynamic_cast<const MeshVectorField *>( _vectorfield );
	    //if( mvf )
	    //build_meshvectorfield_plot();
	    //else
	    build_vectorfield_plot();
	}
    }

    _first = false;
    _oview = _view;
    _olevel = _level;

    if( _field_type != FIELD_NONE && _enabled )
	Colormap::plot( cairo, cm, range );
}


void FieldGraph::plot_sample( cairo_t *cairo, double x, double y, double width, double height )
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


