#include "meshgraph.hpp"
#include "verbose.hpp"


MeshGraph::MeshGraph( const Geometry &g )
    : _g(g)
{

}


MeshGraph::~MeshGraph()
{

}


void MeshGraph::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] )
{
    if( verbose_output )
        std::cout << "  Plotting mesh\n";

    cairo_set_source_rgb( cairo, 0.5, 0.5, 0.5 );
    cairo_set_line_width( cairo, 1.0 );

    int irange[4];
    irange[0] = (int)floor((range[0]-_g.origo(_vb[0])) / _g.h());
    irange[1] = (int)floor((range[1]-_g.origo(_vb[1])) / _g.h());
    irange[2] = (int)ceil((range[2]-_g.origo(_vb[0])) / _g.h());
    irange[3] = (int)ceil((range[3]-_g.origo(_vb[1])) / _g.h());
    
    double xin[2];
    double xout[2];

    for( int i = irange[0]; i <= irange[2]; i++ ) {
	
	//std::cout << "mesh line at x=" << _g.origo(_vb[0]) + _g.h() * i << "\n";
	
	xin[0] = _g.origo(_vb[0]) + _g.h() * i;
	xin[1] = range[1];
	cm->transform( xout, xin );
	cairo_move_to( cairo, xout[0], xout[1] );

	xin[0] = _g.origo(_vb[0]) + _g.h() * i;
	xin[1] = range[3];
	cm->transform( xout, xin );
	cairo_line_to( cairo, xout[0], xout[1] );
    }

    for( int j = irange[1]; j <= irange[3]; j++ ) {
	
	//std::cout << "mesh line at y=" << _g.origo(_vb[1]) + _g.h() * j << "\n";
	
	xin[0] = range[0];
	xin[1] = _g.origo(_vb[1]) + _g.h() * j;
	cm->transform( xout, xin );
	cairo_move_to( cairo, xout[0], xout[1] );

	xin[0] = range[2];
	xin[1] = _g.origo(_vb[1]) + _g.h() * j;
	cm->transform( xout, xin );
	cairo_line_to( cairo, xout[0], xout[1] );
    }

    cairo_stroke( cairo );
}


void MeshGraph::get_bbox( double bbox[4] )
{
    bbox[0] = _g.origo( _vb[0] );
    bbox[1] = _g.origo( _vb[1] );
    bbox[2] = _g.max( _vb[0] );
    bbox[3] = _g.max( _vb[1] );
}

