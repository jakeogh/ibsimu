/*! \file solver2d_smooth.cpp 
 *  \brief Test solver with a 2d surfaces for edge smoothing
 *
 *  \test Test solver with a 2d surfaces for edge smoothing
 */


#include <iostream>
#include <iomanip>
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"

#include "geomplotter.hpp"
#include "frame.hpp"


using namespace std;


bool solid1( double x, double y, double z )
{
    return( x <= 0.025 || y <= 0.025 );
}


bool solid2( double x, double y, double z )
{
    return( x >= 0.045 && x <= 0.055 && y >= 0.055 );
}


void test( int argc, char **argv )
{
    Geometry geom( MODE_2D, Int3D(11,11,1), Vec3D(0,0,0), 0.01 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 3, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, 10.0) );
    geom.build_mesh();

    EpotProblem p;
    //p.enable_smooth_solids( false );
    p.construct( geom );
    //geom.debug_print();
    //p.debug_print();

    ScalarField epot( geom );
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );
    
    GeomPlotter geomplotter( &geom );
    geomplotter.set_mesh( true );
    geomplotter.set_epot( &epot );
    geomplotter.plot_png( "solver2d_smooth.png" );

    /*
    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &g );
    plotter.set_scharge( &scharge );
    plotter.set_epot( &epot );
    plotter.new_geometry_plot_window();
    plotter.run();
    */
}


