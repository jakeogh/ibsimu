/*! \file stl.cpp 
 *  \brief Test STL input.
 *
 *  \test Test STL input.
 */


#include <fstream>
#include <iomanip>
#include "stl_solid.hpp"
#include "stlfile.hpp"
#include "geometry.hpp"
#include "gtkplotter.hpp"
#include "epot_problem.hpp"
#include "error.hpp"
#include "ibsimu.hpp"


using namespace std;


void test( int argc, char **argv )
{
    //Geometry geom( MODE_3D, Int3D(21,21,21), Vec3D(-1.0e-3,-1.0e-3,-1.0e-3), 1e-3 );
    //STLSolid *s1 = new STLSolid( "stl_easy.stl" );

    //Geometry geom( MODE_3D, Int3D(31,71,21), Vec3D(-4e-3,-7e-3,-1e-3), 1e-4 );
    Geometry geom( MODE_3D, Int3D(61,141,41), Vec3D(-4e-3,-7e-3,-1e-3), 5e-5 );
    STLSolid *s1 = new STLSolid( "stl_bin.stl" );
    s1->debug_print( cout );

    geom.set_solid( 7, s1 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 5, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 6, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    geom.build_mesh();



    ScalarField epot( geom );
    ScalarField scharge( geom );

    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.new_geometry_plot_window();
    plotter.run();
}
