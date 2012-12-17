/*! \file surface_inside.cpp
 *  \brief Test surface_inside().
 *
 *  \test Test surface_inside().
 */


#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "geometry.hpp"
#include "func_solid.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "gtkplotter.hpp"


using namespace std;


bool solid1( double x, double y, double z )
{
    return( x < 0.00015 );
}


void test( int argc, char **argv )
{
    double h = 0.001;
    Geometry geom( MODE_3D, Int3D(101,101,101), 
		   Vec3D(-0.05,-0.05,-0.05), h );

    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 3, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 5, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 6, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET, 0.0)  );
    geom.build_mesh();
    geom.build_surface();

    /*
    double x = h/10.0;
    Vec3D v( x, h/10.0, h/10.0 );
    ibsimu.message(1) << x << " " 
		      << geom.inside( v ) << " "
		      << geom.surface_inside( v ) << "\n";
    */
    
    for( double x = -0.002; x < 0.002; x += 0.0001 ) {
	Vec3D v( x, h/10.0, h/10.0 );
	ibsimu.message(1) << x << " " 
			  << geom.inside( v ) << " "
			  << geom.surface_inside( v ) << "\n";
    }

    if( true ) {
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.new_geometry_plot_window();
	plotter.run();
    }
}

