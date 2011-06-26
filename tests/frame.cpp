/*! \file frame.cpp 
 *  \brief Test plot frame.
 *
 *  \test Test plot frame.
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cairo.h>
#include "frame.hpp"
#include "xygraph.hpp"
#include "error.hpp"


using namespace std;


void test( int argc, char **argv )
{
    cairo_surface_t *surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, 640, 480 );
    cairo_t *cairo = cairo_create( surface );

    cairo_set_font_size( cairo, 80 );

    // Make frame
    Frame f;
    f.set_font_size( 15 );
    f.set_geometry( 540, 380, 50, 50 );
    f.set_title( "Title for this plot is very complicated and \nlong "
		 "and it contains \\gamma rke\\epsilon k \\alpha -bet." );
    f.set_axis_label( PLOT_AXIS_X1, "xlabel" );
    f.set_axis_label( PLOT_AXIS_Y1, "ylabel" );
    //f.force_enable_ruler( PLOT_AXIS_Y2, true );
    //f.set_ranges( PLOT_AXIS_Y2, -5.0, 5.0 );

    // Make XYGraph 1
    std::vector<double> xdata;
    std::vector<double> ydata;
    for( size_t a = 0; a < 7; a++ ) {
	xdata.push_back( 100*(-1.0+1.2*a) );
	ydata.push_back( a*a-2.0*a+3.0 );
    }
    XYGraph g( xdata, ydata ); 
    g.set_line_style( XYGRAPH_LINE_SOLID );
    g.set_point_style( XYGRAPH_POINT_CIRCLE, true, 5.0 );
    LegendEntry l( g, "XYGraph 1" );

    // Make XYGraph 2
    std::vector<double> xdata2;
    std::vector<double> ydata2;
    for( size_t a = 0; a < 5; a++ ) {
	xdata2.push_back( -30+182.0*a );
	ydata2.push_back( -2.0*a*a+30.0 );
    }
    XYGraph g2( xdata2, ydata2 ); 
    g2.set_color( Color(0,0,1) );
    g2.set_line_style( XYGRAPH_LINE_SOLID );
    g2.set_point_style( XYGRAPH_POINT_BOX, true, 5.0 );
    LegendEntry l2( g2, "XYGraph boxes" );

    f.add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, &g, &l );
    f.add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, &g2, &l2 );
    f.set_legend_position( LEGEND_POS_TOP_LEFT );
    f.draw( cairo );

    cairo_surface_write_to_png( surface, "frame.png" );
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );
}

