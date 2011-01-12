/*! \file colormap.cpp 
 *  \brief Test colormap.
 *
 *  \test Test colormap.
 */


#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cairo.h>
#if CAIRO_HAS_PS_SURFACE
#include <cairo-ps.h>
#endif
#ifdef CAIRO_HAS_PDF_SURFACE
#include <cairo-pdf.h>
#endif
#include "frame.hpp"
#include "xygraph.hpp"
#include "colormap.hpp"
#include "palette.hpp"
#include "error.hpp"


using namespace std;


void test( void )
{
    cairo_surface_t *surface;
    cairo_t *cairo;

    // Make frame
    Frame f;
    f.set_title( "Title for this plot\nNewline\nAnd another." );
    f.set_axis_label( PLOT_AXIS_X1, "xlabel" );
    f.set_axis_label( PLOT_AXIS_Y1, "ylabel" );
    f.set_ranges( PLOT_AXIS_X1, 0.0, 1.0 );
    f.set_ranges( PLOT_AXIS_Y1, 0.0, 1.0 );

    // Make palette
    std::vector<Palette::Entry> entries;
    entries.push_back( Palette::Entry( Color(1,1,1), 0 ) );
    entries.push_back( Palette::Entry( Color(1,1,0), 1 ) );
    entries.push_back( Palette::Entry( Color(1,0,0), 2 ) );
    entries.push_back( Palette::Entry( Color(0,0,0), 3 ) );
    Palette palette( entries );

    // Make colormap
    double datarange[4] = {0.2, 0.2, 0.8, 0.8};
    std::vector<double> data;
    for( size_t j = 0; j < 15; j++ ) {
	for( size_t i = 0; i < 15; i++ ) {
	    data.push_back( floor(10.0*rand()/(RAND_MAX+1.0)) );
	}
    }
    Colormap colormap( datarange, 15, 15, data );
    colormap.set_palette( palette );
    f.add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, &colormap );

    //
    // PNG
    //

    // Closest interpolation
    colormap.set_interpolation( INTERPOLATION_CLOSEST );
    surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, 640, 480 );
    cairo = cairo_create( surface );
    f.draw( cairo );
    cairo_surface_write_to_png( surface, "colormap_closest.png" );
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );

    // Bilinear interpolation
    colormap.set_interpolation( INTERPOLATION_BILINEAR );
    surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, 640, 480 );
    cairo = cairo_create( surface );
    f.draw( cairo );
    cairo_surface_write_to_png( surface, "colormap_bilinear.png" );
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );

    // Bicubic interpolation
    colormap.set_interpolation( INTERPOLATION_BICUBIC );
    surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, 640, 480 );
    cairo = cairo_create( surface );
    f.draw( cairo );
    cairo_surface_write_to_png( surface, "colormap_bicubic.png" );
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );

    //
    // PS
    //
#ifdef CAIRO_HAS_PS_SURFACE

    // Closest interpolation
    colormap.set_interpolation( INTERPOLATION_CLOSEST );
    surface = cairo_ps_surface_create( "colormap_closest.ps", 640, 480 );
    if( cairo_surface_status( surface ) != CAIRO_STATUS_SUCCESS )
	throw( Error( ERROR_LOCATION, "error creating postscript surface" ) );
#if CAIRO_VERSION >= CAIRO_VERSION_ENCODE(1, 6, 0)
    cairo_ps_surface_set_eps( surface, true );
#endif
    cairo = cairo_create( surface );
    f.draw( cairo );
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );

    // Bilinear interpolation
    colormap.set_interpolation( INTERPOLATION_BILINEAR );
    surface = cairo_ps_surface_create( "colormap_bilinear.ps", 640, 480 );
    if( cairo_surface_status( surface ) != CAIRO_STATUS_SUCCESS )
	throw( Error( ERROR_LOCATION, "error creating postscript surface" ) );
#if CAIRO_VERSION >= CAIRO_VERSION_ENCODE(1, 6, 0)
    cairo_ps_surface_set_eps( surface, true );
#endif
    cairo = cairo_create( surface );
    f.draw( cairo );
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );

    // Bicubic interpolation
    colormap.set_interpolation( INTERPOLATION_BICUBIC );
    surface = cairo_ps_surface_create( "colormap_bicubic.ps", 640, 480 );
    if( cairo_surface_status( surface ) != CAIRO_STATUS_SUCCESS )
	throw( Error( ERROR_LOCATION, "error creating postscript surface" ) );
#if CAIRO_VERSION >= CAIRO_VERSION_ENCODE(1, 6, 0)
    cairo_ps_surface_set_eps( surface, true );
#endif
    cairo = cairo_create( surface );
    f.draw( cairo );
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );

#endif

    //
    // PDF
    //
#ifdef CAIRO_HAS_PDF_SURFACE

    // Closest interpolation
    colormap.set_interpolation( INTERPOLATION_CLOSEST );
    surface = cairo_pdf_surface_create( "colormap_closest.pdf", 640, 480 );
    if( cairo_surface_status( surface ) != CAIRO_STATUS_SUCCESS )
	throw( Error( ERROR_LOCATION, "error creating postscript surface" ) );
    cairo = cairo_create( surface );
    f.draw( cairo );
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );

    // Bilinear interpolation
    colormap.set_interpolation( INTERPOLATION_BILINEAR );
    surface = cairo_pdf_surface_create( "colormap_bilinear.pdf", 640, 480 );
    if( cairo_surface_status( surface ) != CAIRO_STATUS_SUCCESS )
	throw( Error( ERROR_LOCATION, "error creating postscript surface" ) );
    cairo = cairo_create( surface );
    f.draw( cairo );
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );

    // Bicubic interpolation
    colormap.set_interpolation( INTERPOLATION_BICUBIC );
    surface = cairo_pdf_surface_create( "colormap_bicubic.pdf", 640, 480 );
    if( cairo_surface_status( surface ) != CAIRO_STATUS_SUCCESS )
	throw( Error( ERROR_LOCATION, "error creating postscript surface" ) );
    cairo = cairo_create( surface );
    f.draw( cairo );
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );

#endif
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



