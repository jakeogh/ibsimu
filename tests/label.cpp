/*! \file label.cpp 
 *  \test Test label.
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cairo.h>
#include <cmath>
#include "label.hpp"
#include "error.hpp"


using namespace std;


void draw_extents( cairo_t *cairo, double x, double y, 
		   const cairo_text_extents_t *extents )
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


void draw_bbox( cairo_t *cairo, const double bbox[4] )
{
    cairo_set_source_rgb( cairo, 1, 0, 0 );
    cairo_set_line_width( cairo, 1 );
    cairo_rectangle( cairo, bbox[0], bbox[1], 
		     bbox[2]-bbox[0], bbox[3]-bbox[1] );
    cairo_stroke( cairo );
}

void test( void )
{
    cairo_surface_t *surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32,
							   640, 480 );
    cairo_t *cairo = cairo_create( surface );

    cairo_set_source_rgb( cairo, 1, 0, 0 );
    cairo_set_line_width( cairo, 1 );

    cairo_arc( cairo, 100, 100, 2, 0, 2.0*M_PI );
    cairo_fill( cairo );

    cairo_arc( cairo, 200, 100, 2, 0, 2.0*M_PI );
    cairo_fill( cairo );

    cairo_arc( cairo, 300, 100, 2, 0, 2.0*M_PI );
    cairo_fill( cairo );

    cairo_arc( cairo, 100, 200, 2, 0, 2.0*M_PI );
    cairo_fill( cairo );

    cairo_arc( cairo, 200, 200, 2, 0, 2.0*M_PI );
    cairo_fill( cairo );

    cairo_arc( cairo, 300, 200, 2, 0, 2.0*M_PI );
    cairo_fill( cairo );

    cairo_arc( cairo, 300, 300, 2, 0, 2.0*M_PI );
    cairo_fill( cairo );

    cairo_arc( cairo, 200, 300, 2, 0, 2.0*M_PI );
    cairo_fill( cairo );

    cairo_arc( cairo, 100, 300, 2, 0, 2.0*M_PI );
    cairo_fill( cairo );

    cairo_arc( cairo, 300, 400, 2, 0, 2.0*M_PI );
    cairo_fill( cairo );

    cairo_arc( cairo, 200, 400, 2, 0, 2.0*M_PI );
    cairo_fill( cairo );

    cairo_arc( cairo, 100, 400, 2, 0, 2.0*M_PI );
    cairo_fill( cairo );

    double bbox[4];
    Label label;

    // Misc row
    label.set_text( "yksi" );
    label.set_location( 100, 100 );
    label.set_alignment( 0.5, 1.0 );
    label.set_font_size( 50.0 );
    label.draw( cairo );
    label.get_bbox( cairo, bbox );
    draw_bbox( cairo, bbox );

    label.set_text( "kaksi" );
    label.set_location( 300, 100 );
    //label2.set_rotation( -0.5*M_PI );
    label.set_rotation( 0.5 );
    label.set_alignment( 0.5, 1.0 );
    label.set_font_size( 50.0 );
    label.draw( cairo );
    label.get_bbox( cairo, bbox );
    draw_bbox( cairo, bbox );

    // zeroext, Up
    label.set_rotation( 0.0 );
    label.set_text( "kylme\nkala" );
    label.set_location( 300, 200 );
    label.set_alignment( 0.5, 1.0, true );
    label.set_font_size( 30.0 );
    label.draw( cairo );
    label.get_bbox( cairo, bbox );
    draw_bbox( cairo, bbox );

    label.set_text( "u" );
    label.set_location( 200, 200 );
    label.set_alignment( 0.5, 1.0, true );
    label.set_font_size( 30.0 );
    label.draw( cairo );
    label.get_bbox( cairo, bbox );
    draw_bbox( cairo, bbox );

    label.set_text( "neljäyY" );
    label.set_location( 100, 200 );
    label.set_alignment( 0.5, 1.0, true );
    label.set_font_size( 30.0 );
    label.draw( cairo );
    label.get_bbox( cairo, bbox );
    draw_bbox( cairo, bbox );

    // zeroext, Center
    label.set_text( "kylme\nkala" );
    label.set_location( 300, 300 );
    label.set_alignment( 0.5, 0.5, true );
    label.set_font_size( 30.0 );
    label.draw( cairo );
    label.get_bbox( cairo, bbox );
    draw_bbox( cairo, bbox );

    label.set_text( "u" );
    label.set_location( 200, 300 );
    label.set_alignment( 0.5, 0.5, true );
    label.set_font_size( 30.0 );
    label.draw( cairo );
    label.get_bbox( cairo, bbox );
    draw_bbox( cairo, bbox );

    label.set_text( "neljäyY" );
    label.set_location( 100, 300 );
    label.set_alignment( 0.5, 0.5, true );
    label.set_font_size( 30.0 );
    label.draw( cairo );
    label.get_bbox( cairo, bbox );
    draw_bbox( cairo, bbox );

    // zeroext, Down
    label.set_text( "kylme\nkala" );
    label.set_location( 300, 400 );
    label.set_alignment( 0.5, 0.0, true );
    label.set_font_size( 30.0 );
    label.draw( cairo );
    label.get_bbox( cairo, bbox );
    draw_bbox( cairo, bbox );

    label.set_text( "u" );
    label.set_location( 200, 400 );
    label.set_alignment( 0.5, 0.0, true );
    label.set_font_size( 30.0 );
    label.draw( cairo );
    label.get_bbox( cairo, bbox );
    draw_bbox( cairo, bbox );

    label.set_text( "neljäyY" );
    label.set_location( 100, 400 );
    label.set_alignment( 0.5, 0.0, true );
    label.set_font_size( 30.0 );
    label.draw( cairo );
    label.get_bbox( cairo, bbox );
    draw_bbox( cairo, bbox );


    /*
    //Label label( "Toj \u03B1" );
    Label label( "\\alpha=0.0323\n\\beta=0.43331" );
    label.set_location( 200, 200 );
    label.set_alignment( 0.2, 0.75 );
    //label.set_rotation( 3.14159265/6.0 );
    label.set_font_size( 75.0 );
    //label.set_font_weight( CAIRO_FONT_WEIGHT_BOLD );
    //label.set_font_family( "Helvetica" );
    label.draw( cairo );
    cairo_text_extents_t extents;
    label.get_extents( cairo, &extents );
    double bbox[4];
    label.get_bbox( cairo, bbox );
    //draw_extents( cairo, 200, 200, &extents );
    draw_bbox( cairo, bbox );
    */

    cairo_surface_write_to_png( surface, "label.png" );
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



