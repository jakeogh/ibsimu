#include "fieldgraph.hpp"
#include "verbose.hpp"


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
    if( verbose_output )
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













