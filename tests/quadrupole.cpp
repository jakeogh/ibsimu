/*! \file quadrupole.cpp 
 *  \brief Test with a quadrupole field in 2D.
 *
 *  \test Test with a quadrupole field in 2D.
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "geomplotter.hpp"
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


void test( int argc, char **argv )
{
    Geometry geom( MODE_2D, Int3D(11,11,1), Vec3D(0,0,0), 0.01 );

    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    Solid *s3 = new FuncSolid( solid3 );
    geom.set_solid( 9, s3 );
    Solid *s4 = new FuncSolid( solid4 );
    geom.set_solid( 10, s4 );
    geom.set_boundary( 1,  Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 2,  Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 3,  Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4,  Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7,  Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 8,  Bound(BOUND_DIRICHLET, 10.0) );
    geom.set_boundary( 9,  Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 10, Bound(BOUND_DIRICHLET, 10.0) );
    geom.build_mesh();

    EpotProblem p;
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    GeomPlotter geomplotter( &geom );
    geomplotter.set_epot( &epot );
    geomplotter.set_mesh( true );
    geomplotter.plot_png( "quadrupole.png" );
}


