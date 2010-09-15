/*! \file frame_extend_range.cpp 
 *  \brief Test plot frame.
 *
 *  \test Test plot frame.
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cairo.h>
#include "frame.hpp"
#include "xygraph.hpp"
#include "error.hpp"


using namespace std;


void test( void )
{
    cairo_surface_t *surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, 640, 480 );
    cairo_t *cairo = cairo_create( surface );

    // Make frame
    Frame f;
    f.set_geometry( 540, 380, 50, 50 );
    f.set_title( "Title for this plot\nNewline\nAnd another." );
    f.set_axis_label( PLOT_AXIS_X1, "xlabel" );
    f.set_axis_label( PLOT_AXIS_Y1, "ylabel" );
    f.set_fixed_aspect( PLOT_FIXED_ASPECT_EXTEND_RANGE );
    //f.set_fixed_aspect( PLOT_FIXED_ASPECT_INCREASE_MARGIN );
    //f.force_enable_ruler( PLOT_AXIS_Y2, true );
    //f.set_ranges( PLOT_AXIS_Y2, -5.0, 5.0 );

    // Make XYGraph
    std::vector<double> xdata;
    std::vector<double> ydata;
    for( size_t a = 0; a < 21; a++ ) {
	xdata.push_back( 10.0+3.0*sin(2.0*M_PI*a/20.0) );
	ydata.push_back( 20.0+3.0*cos(2.0*M_PI*a/20.0) );
    }
    XYGraph g( xdata, ydata );
    g.set_line_style( XYGRAPH_LINE_DISABLE );
    g.set_point_style( XYGRAPH_POINT_CIRCLE, false, 3.0 );
    f.add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, &g );
    f.draw( cairo );

    cairo_surface_write_to_png( surface, "frame_extend_range.png" );
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



