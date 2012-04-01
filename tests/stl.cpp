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
#include "error.hpp"
#include "ibsimu.hpp"


using namespace std;


void test( int argc, char **argv )
{
    //Geometry geom( MODE_3D, Int3D(41,41,41), Vec3D(0,0,0), 1e-3 );
    //STLSolid *s1 = new STLSolid( "box.stl" );
    
    //Geometry geom( MODE_3D, Int3D(8,8,8), Vec3D(-0.2,-0.2,-0.2), 0.2 );
    //STLSolid *s1 = new STLSolid( "tetra.stl" );

    Geometry geom( MODE_3D, Int3D(41,41,23), Vec3D(10e-3,0e-3,-5e-3), 5e-4 );
    STLSolid *s1 = new STLSolid( "cylinder_low_res_fillet.stl" );

    //Geometry geom( MODE_3D, Int3D(41,41,23), Vec3D(10e-3,0e-3,-5e-3), 5e-4 );
    //STLSolid *s1 = new STLSolid( "cylinder_low_res.stl" );

    geom.set_solid( 7, s1 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 5, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 6, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    geom.build_mesh();

    EpotField epot( geom );

    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.new_geometry_plot_window();
    plotter.run();
}
