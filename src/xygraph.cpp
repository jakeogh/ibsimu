#include <limits>
#include <iostream>
#include <cmath>
#include "xygraph.hpp"


XYGraph::XYGraph()
    : _linewidth(1.0), _color(Color(1,0,0)), _linestyle(XYGRAPH_LINE_DISABLE), _pointstyle(XYGRAPH_POINT_CIRCLE),
      _point_filled(true), _point_scale(3.0)
{
}

XYGraph::XYGraph( const std::vector<double> &xdata, const std::vector<double> &ydata )
    : _linewidth(1.0), _color(Color(1,0,0)), _linestyle(XYGRAPH_LINE_DISABLE), _pointstyle(XYGRAPH_POINT_CIRCLE),
      _point_filled(true), _point_scale(3.0), _xdata(xdata), _ydata(ydata)    
{
}


void XYGraph::plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] ) 
{
    size_t N = _xdata.size() < _ydata.size() ? _xdata.size() : _ydata.size();

    cairo_set_line_width( cairo, _linewidth );
    cairo_set_source_rgba( cairo, _color[0], _color[1], _color[2], _color[3] );

    if( _linestyle == XYGRAPH_LINE_SOLID ) {
	bool cont = false;
	for( size_t a = 0; a < N; a++ ) {
	    double x = _xdata[a];
	    double y = _ydata[a];
	    if( isnan(x) || isnan(y) ) {
		cont = false;
		continue;
	    }
	    cm->transform( x, y );
	    if( cont )
		cairo_line_to( cairo, x, y );
	    else
		cairo_move_to( cairo, x, y );
	    cont = true;
	}
	cairo_stroke( cairo );
    }
    if( _pointstyle != XYGRAPH_POINT_DISABLE ) {
	for( size_t a = 0; a < N; a++ ) {
	    double x = _xdata[a];
	    double y = _ydata[a];
	    if( isnan(x) || isnan(y) )
		continue;

	    cm->transform( x, y );
	    cairo_save( cairo );
	    cairo_translate( cairo, x, y );
	    cairo_scale( cairo, _point_scale, _point_scale );

	    if( _pointstyle == XYGRAPH_POINT_CIRCLE ) {
		cairo_move_to( cairo, 1.0, 0.0 );
		cairo_arc( cairo, 0.0, 0.0, 1.0, 0.0, 2.0*M_PI );
	    }

	    cairo_restore( cairo );
	    if( _point_filled )
		cairo_fill( cairo );
	    else
		cairo_stroke( cairo );
	}
    }
}


void XYGraph::get_bbox( double bbox[4] )
{
    size_t N = _xdata.size() < _ydata.size() ? _xdata.size() : _ydata.size();

    bbox[0] = std::numeric_limits<double>::infinity();
    bbox[1] = std::numeric_limits<double>::infinity();
    bbox[2] = -std::numeric_limits<double>::infinity();
    bbox[3] = -std::numeric_limits<double>::infinity();
    for( size_t a = 0; a < N; a++ ) {
	if( _xdata[a] < bbox[0] )
	    bbox[0] = _xdata[a];
	if( _xdata[a] > bbox[2] )
	    bbox[2] = _xdata[a];
	if( _ydata[a] < bbox[1] )
	    bbox[1] = _ydata[a];
	if( _ydata[a] > bbox[3] )
	    bbox[3] = _ydata[a];
    }
}


void XYGraph::set_data( const std::vector<double> &xdata, 
			const std::vector<double> &ydata )
{
    _xdata = xdata;
    _ydata = ydata;
}


void XYGraph::set_color( const Color &color )
{
    _color = color;
}


void XYGraph::set_line_width( double linewidth )
{
    _linewidth = linewidth;
}


void XYGraph::set_line_style( line_style_e linestyle )
{
    _linestyle  = linestyle;
}


void XYGraph::set_point_style( point_style_e pointstyle, bool filled, double scale )
{
    _pointstyle = pointstyle;
    _point_filled = filled;
    _point_scale = scale;
}













