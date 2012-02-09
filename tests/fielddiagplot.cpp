/*! \file fielddiagplot.cpp 
 *  \brief Test field diagostic plotter.
 *
 *  \test Test field diagostic plotter.
 *
 */


#include <sstream>
#include <fstream>
#include <iomanip>
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "particledatabase.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "vectorfield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "particlediagplotter.hpp"
#include "fielddiagplotter.hpp"


using namespace std;


bool solid1( double x, double y, double z )
{
    return( x <= 0.00187 && y >= 0.00054 && y >= 2.28*x - 0.0010 &&
	    (x >= 0.00054 || y >= 0.0015) );
}


bool solid2( double x, double y, double z )
{
    return( x >= 0.0095 && y >= 0.0023333 && y >= 0.01283 - x );
}


void test( int argc, char **argv )
{
    Geometry geom( MODE_2D, Int3D(241,141,1), Vec3D(0,0,0), 0.00005 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0 ) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, -8.0e3) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0)  );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, -8.0e3) );
    geom.build_mesh();

    EpotProblem p;
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    /*
    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_scharge( &scharge );
    plotter.set_epot( &epot );
    plotter.new_geometry_plot_window();
    plotter.run();
    */

    FieldDiagPlotter fplotter( geom );
    fplotter.set_scharge( &scharge );
    fplotter.set_epot( &epot );
    fplotter.set_coordinates( 100, Vec3D(0,0,0), Vec3D(12e-3,0,0) );

    field_diag_type_e diag[2] = {FIELDD_DIAG_EPOT, FIELDD_DIAG_NONE};
    field_loc_type_e loc[2] = {FIELDD_LOC_X, FIELDD_LOC_NONE};
    fplotter.set_diagnostic( diag, loc );
    fplotter.plot_png( "fielddiagplot.png" );

    field_diag_type_e diag2[2] = {FIELDD_DIAG_SCHARGE, FIELDD_DIAG_NONE};
    field_loc_type_e loc2[2] = {FIELDD_LOC_X, FIELDD_LOC_NONE};
    fplotter.set_diagnostic( diag2, loc2 );
    fplotter.plot_png( "fielddiagplot2.png" );
}


