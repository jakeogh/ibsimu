/*! \file scharge_out_of_range.cpp 
 *  \brief Test for scharge out of range error
 *
 *  \test Test for scharge out of range error
 *
 */


#include <iostream>
#include <iomanip>
#include "geometry.hpp"
#include "particledatabase.hpp"
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "epot_efield.hpp"
#include "particles.hpp"
#include "error.hpp"
#include "verbose.hpp"

#include "gtkplotter.hpp"

using namespace std;


void test( int *argc, char ***argv )
{
    verbose_output = 1;

    Geometry g( MODE_2D, Int3D(11,11,1), Vec3D(0,0,0), 0.01 );
    g.set_boundary( 1, Bound(BOUND_DIRICHLET, -1000.0) );
    g.set_boundary( 2, Bound(BOUND_DIRICHLET,     0.0) );
    g.set_boundary( 3, Bound(BOUND_NEUMANN,       0.0) );
    g.set_boundary( 4, Bound(BOUND_NEUMANN,       0.0) );
    g.build_mesh();

    EpotProblem p;
    p.construct( g );

    ScalarField epot( g );
    VectorField bfield;
    ScalarField scharge( g );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    EpotEfield efield( g, epot );
    efield_extrpl_e extrpl[6] = { EFIELD_MIRROR, EFIELD_MIRROR,
				  EFIELD_MIRROR, EFIELD_MIRROR,
				  EFIELD_MIRROR, EFIELD_MIRROR };
    efield.set_extrapolation( extrpl );

    ParticleDataBase2D pdb;
    pdb.set_thread_count( 1 );
    //bool pmirror[6] = { false, false, false, false, false, false };
    bool pmirror[6] = { false, true, false, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0, 0.1, 1e4, 0.0, 1e5 ) );
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0, 0.05, 5e5, 0.05, 1e5 ) );
    pdb.iterate_trajectories( scharge, efield, bfield, g );
    pdb.debug_print();

    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &g );
    plotter.set_scharge( &scharge );
    plotter.set_epot( &epot );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();

    /*
    GeomPlotter geomplotter;
    geomplotter.set_geometry( g );
    geomplotter.set_epot( epot );
    geomplotter.set_particledatabase( pdb );
    pngplot( &geomplotter, "scharge_out_of_range.png" );
    */
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

