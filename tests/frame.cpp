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


void test( void )
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

    // Make XYGraph
    std::vector<double> xdata;
    std::vector<double> ydata;
    for( size_t a = 0; a < 7; a++ ) {
	xdata.push_back( 100*(-1.0+1.2*a) );
	ydata.push_back( a*a-2.0*a+3.0 );
    }
    XYGraph g( xdata, ydata ); 
    g.set_line_style( XYGRAPH_LINE_SOLID );
    g.set_point_style( XYGRAPH_POINT_CIRCLE, true, 5.0 );
    f.add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, &g );
    f.draw( cairo );

    cairo_surface_write_to_png( surface, "frame.png" );
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );
}


int main( void )
{
    try {
	test();
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

