/*! \file fieldgraph.cpp 
 *  \brief Test fieldgraph
 *
 *  \test Test fieldgraph
 */

#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "geomplotter.hpp"


bool solid1( double x, double y, double z )
{
    return( x < -0.05 && x > -0.06 && y < 0.08 && y > -0.08 );
}


bool solid2( double x, double y, double z )
{
    return( x < 0.06 && x > 0.05 && y < 0.08 && y > -0.08 );
}


void test( int argc, char **argv )
{
    Geometry geom( MODE_2D, Int3D(121,121,1), Vec3D(-0.12,-0.12,0), 0.002 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary( 3, Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET, -10.0) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET,  10.0) );
    geom.build_mesh();

    EpotProblem p;
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );
    EpotEfield efield( geom, epot );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    GeomPlotter gplotter( &geom );
    gplotter.set_epot( &epot );
    gplotter.set_efield( &efield );
    gplotter.set_size( 800, 600 );
    gplotter.set_font_size( 16 );

    gplotter.set_fieldgraph_plot( FIELD_EFIELD );
    gplotter.plot_png( "fieldgraph_efield.png" );

    gplotter.set_fieldgraph_plot( FIELD_EPOT );
    gplotter.plot_png( "fieldgraph_epot.png" );

}
