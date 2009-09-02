/*! \file gtkplotter.cpp 
 *  \brief Test gtkplotter.
 *
 *  \test Test gtkplotter.
 */


#include <iostream>
#include <iomanip>
#include "func_solid.hpp"
#include "gtkplotter.hpp"
#include "error.hpp"


using namespace std;


bool s1( double x, double y, double z )
{
    return( x*x + y*y < 0.1*0.1 );
}

void test( int *argc, char ***argv )
{
    Geometry geom( MODE_2D, Int3D(20,20,1), Vec3D(0,0,0), 0.01 );
    Solid *s = new FuncSolid( s1 );
    geom.set_solid( 7, s );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET, 20.0) );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET, 20.0) );
    geom.build_mesh();

    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.new_geometry_plot_window();
    plotter.run();
}


int main( int argc, char **argv )
{
    try {
	test( &argc, &argv );
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

