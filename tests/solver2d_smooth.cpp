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

#include "gtkplotter.hpp"
#include "frame.hpp"
#include "solidplot.hpp"
#include "particleplot.hpp"
#include "eqpotplot.hpp"

using namespace std;


bool solid1( double x, double y, double z )
{
    return( x <= 0.025 || y <= 0.025 );
}


bool solid2( double x, double y, double z )
{
    return( x >= 0.045 && x <= 0.055 && y >= 0.055 );
}


void test( int *argc, char ***argv )
{
    verbose_output = 1;

    Geometry g( MODE_2D, Int3D(11,11,1), Vec3D(0,0,0), 0.01 );
    Solid *s1 = new FuncSolid( solid1 );
    g.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    g.set_solid( 8, s2 );
    g.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 2, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 3, Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 8, Bound(BOUND_DIRICHLET, 10.0) );
    g.build_mesh();

    EpotProblem p;
    //p.enable_smooth_solids( false );
    p.construct( g );
    g.debug_print();
    p.debug_print();

    ScalarField epot( g );
    ScalarField scharge( g );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );
    
    /*
    GeomPlotter geomplotter;
    geomplotter.set_size( 1024, 768 );
    //geomplotter.set_range( 0.09, 0.12000001, 0.02, 0.04 );
    geomplotter.set_meshlines( true );
    geomplotter.set_geometry( g );
    geomplotter.set_epot( epot );
    pngplot( &geomplotter, "solver2d_smooth.png" );
    */

    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &g );
    plotter.set_scharge( &scharge );
    plotter.set_epot( &epot );
    //plotter.set_particledatabase( &pdb );
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

