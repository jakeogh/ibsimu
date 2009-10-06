/*! \file fonts.cpp 
 *  \test Test fonts.
 */


#include <iostream>
#include <iomanip>
#include <cairo.h>
#include "fonts.hpp"
#include "error.hpp"


using namespace std;


void draw_extents( cairo_t *cairo, double x, double y, const cairo_text_extents_t *extents )
{
    cairo_set_source_rgb( cairo, 1, 0, 0 );
    cairo_set_line_width( cairo, 1 );

    cairo_move_to( cairo, 
		   x+extents->x_bearing, 
		   y );
    cairo_line_to( cairo,
		   x+extents->x_bearing+extents->width, 
		   y );
    cairo_move_to( cairo,
		   x, 
		   y+extents->y_bearing );
    cairo_line_to( cairo,
		   x, 
		   y+extents->y_bearing+extents->height );
    cairo_move_to( cairo, 
		   x+extents->x_bearing, 
		   y+extents->y_bearing+extents->height );
    cairo_line_to( cairo, 
		   x+extents->x_bearing+extents->width, 
		   y+extents->y_bearing+extents->height );
    cairo_line_to( cairo, 
		   x+extents->x_bearing+extents->width, 
		   y+extents->y_bearing );
    cairo_line_to( cairo, 
		   x+extents->x_bearing, 
		   y+extents->y_bearing );
    cairo_line_to( cairo, 
		   x+extents->x_bearing, 
		   y+extents->y_bearing+extents->height );
  
    cairo_stroke( cairo );
}


void test( void )
{
    cairo_surface_t *surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, 640, 480 );
    cairo_t *cairo = cairo_create( surface );
    cairo_set_font_size( cairo, 20 );

    cairo_text_extents_t extents;
    cairo_matrix_t matrix;

    cairo_set_font_size( cairo, 100.0 );
    cairo_get_font_matrix( cairo, &matrix );
    cairo_matrix_rotate( &matrix, -0.1 );
    cairo_set_font_matrix( cairo, &matrix );

    double x = 100;
    double y = 200;
    fontlib.draw_text( cairo, "Toj \u03B1 \u2630", x, y );
    fontlib.text_extents( cairo, "Toj \u03B1 \u2630", &extents );
    draw_extents( cairo, 100, 200, &extents );

    fontlib.draw_text( cairo, "Kala", x, y );

    cairo_surface_write_to_png( surface, "fonts.png" );
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

