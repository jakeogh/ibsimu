/*! \file transformation.cpp 
 *  \brief Test Transformation class.
 *
 *  \test Test Transformation class.
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cairo.h>
#include <math.h>
#include "frame.hpp"
#include "xygraph.hpp"
#include "transformation.hpp"
#include "error.hpp"


using namespace std;


class Bitmap : public Graph {
    
    const Transformation &_T;

public:

    Bitmap( const Transformation &T ) : _T(T) {}

    virtual ~Bitmap() {}

    virtual void plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] ) {

	// Calculate pixel ranges
	double prange[4];
	cm->transform( &prange[0], &range[0] );
	cm->transform( &prange[2], &range[2] );

	// Calculate pixel integer limits of drawn area, 
	// y flipped to have smaller numbers as 0 and 1, 
	// bigger as 2 and 3.
	int plim[4] = { (int)floor(prange[0]+0.5),
			(int)floor(prange[3]+0.5),
			(int)floor(prange[2]+0.5),
			(int)floor(prange[1]+0.5) };

	cairo_surface_t *surface = cairo_get_target( cairo );
	unsigned char *buf = cairo_image_surface_get_data( surface );
	int width  = cairo_image_surface_get_width( surface );
	int height = cairo_image_surface_get_height( surface );
	int stride = cairo_image_surface_get_stride( surface );

	if( plim[0] < 0 )
	    plim[0] = 0;
	if( plim[1] < 0 )
	    plim[1] = 0;
	if( plim[2] >= width )
	    plim[2] = width-1;
	if( plim[3] >= height )
	    plim[3] = height-1;

	for( int i = plim[0]; i <= plim[2]; i++ ) {
	    for( int j = plim[1]; j <= plim[3]; j++ ) {

		double x[2] = { i, j };
		cm->inv_transform( x[0], x[1] );

		Vec3D v( x[0], x[1] );
		v = _T.inv_transform_point( v );

		Color c;
		if( (v[0] > -2.0 && v[0] < 2.0 &&
		     v[1] > -2.0 && v[1] < 2.0) && 
		    !(v[0] > 1.0 && v[1] < 0.0) )
		    c = Color(0,0,1);
		else
		    c = Color(1,1,1);

		buf[j*stride+4*i+0] = (unsigned char)(255*c[2]);  // Blue
		buf[j*stride+4*i+1] = (unsigned char)(255*c[1]);  // Green
		buf[j*stride+4*i+2] = (unsigned char)(255*c[0]);  // Red
		buf[j*stride+4*i+3] = (unsigned char)255;         // Alpha
	    }
	}
    }

    virtual void get_bbox( double bbox[4] ) {
	bbox[0] = 0.0;
	bbox[1] = 0.0;
	bbox[2] = 0.0;
	bbox[3] = 0.0;
    }
};


void plot( const Transformation &T, const std::string &outfile, const std::string &title )
{
    cairo_surface_t *surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, 640, 480 );
    cairo_t *cairo = cairo_create( surface );

    // Make frame
    Frame f;
    f.set_fixed_aspect( PLOT_FIXED_ASPECT_EXTEND_RANGE );
    f.set_font_size( 12 );
    f.set_geometry( 640, 480, 0, 0 );
    f.set_ranges( PLOT_AXIS_X1, -8.0, 8.0 );
    f.set_ranges( PLOT_AXIS_Y1, -8.0, 8.0 );

    // Make XYGraph
    std::vector<double> xdata;
    std::vector<double> ydata;

    xdata.push_back( -2.0 );
    ydata.push_back( -2.0 );

    xdata.push_back(  1.0 );
    ydata.push_back( -2.0 );

    xdata.push_back(  1.0 );
    ydata.push_back(  0.0 );

    xdata.push_back(  2.0 );
    ydata.push_back(  0.0 );

    xdata.push_back(  2.0 );
    ydata.push_back(  2.0 );

    xdata.push_back( -2.0 );
    ydata.push_back(  2.0 );

    xdata.push_back( -2.0 );
    ydata.push_back( -2.0 );

    // Transform xygraph (direct)
    for( size_t a = 0; a < xdata.size(); a++ ) {
	Vec3D v( xdata[a], ydata[a] );
	v = T.transform_point( v );
	xdata[a] = v[0];
	ydata[a] = v[1];
    }

    Bitmap bitmap( T );
    f.add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, &bitmap );

    XYGraph box( xdata, ydata ); 
    box.set_line_style( XYGRAPH_LINE_SOLID );
    box.set_line_width( 3.0 );
    box.set_color( Color(0,0,0) );
    box.set_point_style( XYGRAPH_POINT_DISABLE, true, 1.0 );
    f.add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, &box );

    f.set_title( title );
    f.draw( cairo );

    cairo_surface_write_to_png( surface, outfile.c_str() );
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );
}


void test( int argc, char **argv )
{
    Transformation T;
    plot( T, "transformation1.png", "No transformation" );

    T.translate( Vec3D(5,0,0) );
    plot( T, "transformation2.png", "Translation" );

    T.rotate_z( 45.0*M_PI/180.0 );
    plot( T, "transformation3.png", "Translation + rotation" );

    T.reset();

    T.rotate_z( 45.0*M_PI/180.0 );
    plot( T, "transformation4.png", "Rotation" );

    T.translate( Vec3D(5,0,0) );
    plot( T, "transformation5.png", "Rotation + translation" );
}

