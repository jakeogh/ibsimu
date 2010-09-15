/*! \file ruler_inv.cpp 
 *  \brief Test plot ruler inversed min max ranges.
 *
 *  \test Test plot ruler inversed min max ranges.
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cairo.h>
#include "ruler.hpp"
#include "error.hpp"


using namespace std;


void draw_bbox( cairo_t *cairo, const double bbox[4] )
{
    cairo_set_source_rgb( cairo, 1, 0, 0 );
    cairo_set_line_width( cairo, 1 );
    cairo_rectangle( cairo, bbox[0], bbox[1], bbox[2]-bbox[0], bbox[3]-bbox[1] );
    cairo_stroke( cairo );
}


void test( void )
{
    cairo_surface_t *surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, 640, 480 );
    cairo_t *cairo = cairo_create( surface );
    cairo_set_font_size( cairo, 20 );

    Coordmapper1D cm;
    double bbox[4];

    Ruler xruler( 0 );
    xruler.set_font_size( 20 );
    xruler.set_endpoints( 80, 400, 550, 400 );
    xruler.set_ranges( 1.2, -2.3 );
    xruler.set_ticlen( 5.0, 10.0 );
    xruler.set_axis_label( "xlabel" );
    xruler.draw( cairo, cm );
    xruler.get_bbox( cairo, bbox, cm );
    draw_bbox( cairo, bbox );

    cairo_surface_write_to_png( surface, "ruler_inv.png" );
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



