#include <iostream>
#include <iomanip>
#include "fileplot.hpp"
#include "gtkplot.hpp"
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"


using namespace std;


const double r = 0.04;


bool solid1( double x, double y, double z )
{
    return( x*x+y*y <= r*r );
}

bool solid2( double x, double y, double z )
{
    x -= 0.1;
    return( x*x+y*y <= r*r );
}

bool solid3( double x, double y, double z )
{
    x -= 0.1;
    y -= 0.1;
    return( x*x+y*y <= r*r );
}

bool solid4( double x, double y, double z )
{
    y -= 0.1;
    return( x*x+y*y <= r*r );
}


void test( int *argc, char ***argv )
{
    //verbose_output = 1;

    // size = (0.1,0.1)
    //Geometry g( MODE_2D, Int3D(51,51,1), Vec3D(0,0,0), 0.002 );
    //Geometry g( MODE_2D, Int3D(26,26,1), Vec3D(0,0,0), 0.004 );
    Geometry g( MODE_2D, Int3D(11,11,1), Vec3D(0,0,0), 0.01 );

    Solid *s1 = new FuncSolid( solid1 );
    g.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    g.set_solid( 8, s2 );
    Solid *s3 = new FuncSolid( solid3 );
    g.set_solid( 9, s3 );
    Solid *s4 = new FuncSolid( solid4 );
    g.set_solid( 10, s4 );
    g.set_boundary( 1,  Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 2,  Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 3,  Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 4,  Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 7,  Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 8,  Bound(BOUND_DIRICHLET, 10.0) );
    g.set_boundary( 9,  Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 10, Bound(BOUND_DIRICHLET, 10.0) );
    g.build_mesh();

    EpotProblem p;
    p.construct( g );
    //g.debug_print();
    //p.debug_print();

    ScalarField epot( g );
    ScalarField scharge( g );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    GeomPlotter geomplotter;
    geomplotter.set_geometry( g );
    geomplotter.set_epot( epot );
    geomplotter.set_meshlines( true );
    //gtkplot( &geomplotter, argc, argv );
    pngplot( &geomplotter, "quadrupole.png" );
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

