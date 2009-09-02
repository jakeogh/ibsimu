#include <limits>
#include "geomplot.hpp"


GeomPlot::GeomPlot( Frame *frame, const Geometry *geom )
    : _frame(frame), _geom(geom), _epot(NULL), _scharge(NULL), _pdb(NULL),
      _solidgraph(NULL), _fieldgraph(NULL), _eqpotgraph(NULL), _particlegraph(NULL),
      _meshgraph(NULL), _view(VIEW_XY), _level(0),
      _eqlines_auto(20), _particle_div(10), _scharge_field(false), _qm_discretation(true), 
      _mesh(false), _cache(true)
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
    if( _particlegraph )
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _particlegraph );
    if( _fieldgraph )
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _fieldgraph );
    if( _solidgraph )
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _solidgraph );
    if( _eqpotgraph )
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _eqpotgraph );
    if( _meshgraph )
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _meshgraph );
}


void GeomPlot::set_epot( const ScalarField *epot ) 
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
	_eqpotgraph->set_eqlines_auto( _eqlines_auto );
	_eqpotgraph->set_eqlines_manual( _eqlines_manual );
    }

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


void GeomPlot::set_scharge( const ScalarField *scharge ) 
{
    _scharge = scharge;
    // space charge field graph is disabled by default
}


void GeomPlot::set_scharge_field( bool enable ) 
{
    _scharge_field = enable;

    if( _fieldgraph ) {
	delete _fieldgraph;
	_fieldgraph = NULL;
    }

    if( _scharge && _scharge_field ) {
	_fieldgraph = new FieldGraph( _scharge );
    }

    reset_graphs();
}



void GeomPlot::set_particle_database( const ParticleDataBase *pdb ) 
{
    _pdb = pdb;

    if( _particlegraph ) {
	delete _particlegraph;
	_particlegraph = NULL;
    }

    if( _pdb ) {
	_particlegraph = new ParticleGraph( *_geom, *_pdb );
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
    }

    reset_graphs();
}


void GeomPlot::set_view( view_e view, int level ) 
{
    if( view == VIEW_XY ) {
	_vb[0] = 0;
	_vb[1] = 1;
	_vb[2] = 2;

	// Set axis labels
	_frame->set_axis_label( PLOT_AXIS_X1, "x (m)" );
	if( _geom->geom_mode() == MODE_CYL )
	    _frame->set_axis_label( PLOT_AXIS_Y1, "r (m)" );
	else
	    _frame->set_axis_label( PLOT_AXIS_Y1, "y (m)" );

    } else if( view == VIEW_XZ ) {
	_vb[0] = 0;
	_vb[1] = 2;
	_vb[2] = 1;

	// Set axis labels
	_frame->set_axis_label( PLOT_AXIS_X1, "x (m)" );
	_frame->set_axis_label( PLOT_AXIS_Y1, "z (m)" );

    } else if( view == VIEW_YZ ) {
	_vb[0] = 1;
	_vb[1] = 2;
	_vb[2] = 0;

	// Set axis labels
	if( _geom->geom_mode() == MODE_CYL )
	    _frame->set_axis_label( PLOT_AXIS_X1, "r (m)" );
	else
	    _frame->set_axis_label( PLOT_AXIS_X1, "y (m)" );
	_frame->set_axis_label( PLOT_AXIS_Y1, "z (m)" );

    }

    // Set view and level
    _view = view;
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


