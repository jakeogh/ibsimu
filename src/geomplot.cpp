/*! \file geomplot.cpp
 *  \brief Geometry plotting
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
#include "geomplot.hpp"


GeomPlot::GeomPlot( Frame *frame, const Geometry *geom )
    : _frame(frame), _geom(geom), _epot(NULL), _scharge(NULL), _tdens(NULL), _bfield(NULL),
      _pdb(NULL), _solidgraph(NULL), _fieldgraph(NULL), _eqpotgraph(NULL), _particlegraph(NULL),
      _meshgraph(NULL), _view(VIEW_XY), _level(0),
      _eqlines_auto(20), _particle_div(10), _scharge_field(false), _qm_discretation(true), 
      _mesh(false), _fieldplot_sel(FIELD_NONE), _fieldplot_logscale(false), _cache(true)
{
    if( _geom == NULL )
	throw( Error( ERROR_LOCATION, "geometry undefined" ) );

    // Set frame basic properties
    _frame->set_fixed_aspect( PLOT_FIXED_ASPECT_INCREASE_MARGIN );
    double min = -std::numeric_limits<double>::infinity();
    double max = std::numeric_limits<double>::infinity();
    _frame->set_ranges( PLOT_AXIS_X1, min, max );
    _frame->set_ranges( PLOT_AXIS_Y1, min, max );
    _frame->ruler_autorange_enable( PLOT_AXIS_X1, false, false );
    _frame->ruler_autorange_enable( PLOT_AXIS_Y1, false, false );

    // Add default drawable (solid geometry)
    _solidgraph = new SolidGraph( *_geom );
    _frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _solidgraph );
    if( _geom->geom_mode() == MODE_3D )
	set_view( VIEW_XY, -1 );
    else
	set_view( VIEW_XY, 0 );
}


GeomPlot::~GeomPlot()
{
    if( _solidgraph )
	delete _solidgraph;
    if( _fieldgraph )
	delete _fieldgraph;
    if( _eqpotgraph )
	delete _eqpotgraph;
    if( _particlegraph )
	delete _particlegraph;
    if( _meshgraph )
	delete _meshgraph;
}


void GeomPlot::disable_cache( void )
{
    _cache = false;
    if( _solidgraph )
	_solidgraph->disable_cache();
    if( _eqpotgraph )
	_eqpotgraph->disable_cache();
}


void GeomPlot::reset_graphs()
{
    _frame->clear_graphs();

    // Ensure correct order of graphs
    if( _fieldgraph )
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _fieldgraph );
    if( _particlegraph )
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _particlegraph );
    if( _solidgraph )
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _solidgraph );
    if( _eqpotgraph )
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _eqpotgraph );
    if( _meshgraph )
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _meshgraph );
}


void GeomPlot::set_epot( const MeshScalarField *epot ) 
{
    _epot = epot;

    if( _eqpotgraph ) {
	delete _eqpotgraph;
	_eqpotgraph = NULL;
    }

    if( _epot ) {
	_eqpotgraph = new EqPotGraph( *_epot, *_geom );
	if( !_cache )
	    _eqpotgraph->disable_cache();
	_eqpotgraph->set_view( _view, _level );
	_eqpotgraph->set_eqlines_auto( _eqlines_auto );
	_eqpotgraph->set_eqlines_manual( _eqlines_manual );
    }

    if( _fieldplot_sel == FIELD_EPOT ) {
	// Redo field graph to use new epot
	set_fieldgraph_plot( _fieldplot_sel );
    }
    
    // Reset graphs in every case
    reset_graphs();
}


void GeomPlot::set_eqlines_manual( const std::vector<double> &pot )
{
    _eqlines_manual = pot;
    if( _eqpotgraph )
	_eqpotgraph->set_eqlines_manual( pot );
}


void GeomPlot::set_eqlines_auto( size_t N ) 
{
    _eqlines_auto = N;
    if( _eqpotgraph )
	_eqpotgraph->set_eqlines_auto( N );
}


void GeomPlot::set_trajdens( const MeshScalarField *tdens ) 
{
    _tdens = tdens;
    if( _fieldplot_sel == FIELD_TRAJDENS ) {
	// Redo field graph to use new tdens
	set_fieldgraph_plot( _fieldplot_sel );
    }
}


void GeomPlot::set_scharge( const MeshScalarField *scharge ) 
{
    _scharge = scharge;
    if( _fieldplot_sel == FIELD_SCHARGE ) {
	// Redo field graph to use new scharge
	set_fieldgraph_plot( _fieldplot_sel );
    }
}


void GeomPlot::set_fieldgraph_plot( field_type_e fieldplot )
{
    _fieldplot_sel = fieldplot;

    // Remove old field graph.
    if( _fieldgraph ) {
	delete _fieldgraph;
	_fieldgraph = NULL;
    }

    // Define new field graph.
    switch( _fieldplot_sel ) {
    case FIELD_EPOT:
	if( _epot )
	    _fieldgraph = new FieldGraph( _epot );
	break;
    case FIELD_SCHARGE:
	if( _scharge )
	    _fieldgraph = new FieldGraph( _scharge );
	break;
    case FIELD_TRAJDENS:
	if( _tdens )
	    _fieldgraph = new FieldGraph( _tdens );
	break;
    case FIELD_NONE:
	// No fieldgraph defined
	break;
    default:
	throw( ErrorUnimplemented( ERROR_LOCATION, "Unimplemented field plottingd" ) );
	break;
    }

    if( _fieldgraph ) {
	_fieldgraph->set_view( _view, _level );
	_fieldgraph->set_logscale( _fieldplot_logscale );
    }

    reset_graphs();
}


void GeomPlot::set_fieldgraph_logscale( bool enable )
{
    _fieldplot_logscale = enable;
    if( _fieldgraph ) {
	_fieldgraph->set_logscale( _fieldplot_logscale );
    }
}


void GeomPlot::set_particledatabase( const ParticleDataBase *pdb ) 
{
    _pdb = pdb;

    if( _particlegraph ) {
	delete _particlegraph;
	_particlegraph = NULL;
    }

    if( _pdb ) {
	_particlegraph = new ParticleGraph( *_geom, *_pdb );
	_particlegraph->set_view( _view, _level );
	_particlegraph->set_particle_div( _particle_div );
	_particlegraph->set_qm_discretation( _qm_discretation );
    }

    reset_graphs();
}


void GeomPlot::set_particle_div( size_t particle_div ) 
{
    _particle_div = particle_div;
    if( _particlegraph )
	_particlegraph->set_particle_div( _particle_div );
}


void GeomPlot::set_qm_discretation( bool enable ) 
{
    _qm_discretation = enable;
    if( _particlegraph )
	_particlegraph->set_qm_discretation( _qm_discretation );
}



void GeomPlot::set_mesh( bool enable ) 
{
    _mesh = enable;

    if( _meshgraph ) {
	delete _meshgraph;
	_meshgraph = NULL;
    }

    if( _mesh ) {
	_meshgraph = new MeshGraph( *_geom );
	_meshgraph->set_view( _view, _level );
    }

    reset_graphs();
}


void GeomPlot::set_view( view_e view, int level ) 
{
    //std::cout << "set_view( " << view << ", " << level << " )\n";

    switch( view ) {
    case VIEW_XY:
	_vb[0] = 0;
	_vb[1] = 1;
	_vb[2] = 2;

	// Set axis labels
	_frame->set_axis_label( PLOT_AXIS_X1, "x (m)" );
	if( _geom->geom_mode() == MODE_CYL )
	    _frame->set_axis_label( PLOT_AXIS_Y1, "r (m)" );
	else
	    _frame->set_axis_label( PLOT_AXIS_Y1, "y (m)" );
	break;

    case VIEW_XZ:
	_vb[0] = 0;
	_vb[1] = 2;
	_vb[2] = 1;

	// Set axis labels
	_frame->set_axis_label( PLOT_AXIS_X1, "x (m)" );
	_frame->set_axis_label( PLOT_AXIS_Y1, "z (m)" );
	break;

    case VIEW_YX:
	_vb[0] = 1;
	_vb[1] = 0;
	_vb[2] = 2;

	// Set axis labels
	_frame->set_axis_label( PLOT_AXIS_X1, "y (m)" );
	_frame->set_axis_label( PLOT_AXIS_Y1, "x (m)" );
	break;

    case VIEW_YZ:
	_vb[0] = 1;
	_vb[1] = 2;
	_vb[2] = 0;

	// Set axis labels
	_frame->set_axis_label( PLOT_AXIS_X1, "y (m)" );
	_frame->set_axis_label( PLOT_AXIS_Y1, "z (m)" );
	break;

    case VIEW_ZX:
	_vb[0] = 2;
	_vb[1] = 0;
	_vb[2] = 1;

	// Set axis labels
	_frame->set_axis_label( PLOT_AXIS_X1, "z (m)" );
	_frame->set_axis_label( PLOT_AXIS_Y1, "x (m)" );
	break;

    case VIEW_ZY:
	_vb[0] = 2;
	_vb[1] = 1;
	_vb[2] = 0;

	// Set axis labels
	_frame->set_axis_label( PLOT_AXIS_X1, "z (m)" );
	_frame->set_axis_label( PLOT_AXIS_Y1, "y (m)" );
	break;

    default:
	throw( ErrorUnimplemented( ERROR_LOCATION ) );
	break;
    }

    // Set view
    _view = view;
    
    // Set and check level
    if( level == -1 )
	_level = _geom->size(_vb[2])/2;
    else if( level < 0 )
	_level = 0;
    else if( level >= (int32_t)_geom->size(_vb[2]) )
	_level = _geom->size(_vb[2])-1;
    else
	_level = level;

    if( _solidgraph )
	_solidgraph->set_view( _view, _level );
    if( _eqpotgraph )
	_eqpotgraph->set_view( _view, _level );
    if( _particlegraph )
	_particlegraph->set_view( _view, _level );
    if( _meshgraph )
	_meshgraph->set_view( _view, _level );
    if( _fieldgraph )
	_fieldgraph->set_view( _view, _level );
}







