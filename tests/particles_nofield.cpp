/*! \file particles_nofield.cpp 
 *  \brief Test particle iterator in zero electric field.
 *
 *  \test Test particle iterator in zero electric field.
 */

/* Particle iterator fails using adaptive step size method if there is
 * no fields. The iterator can't estimate the error.
 *
 */

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "particles.hpp"
#include "gtkplotter.hpp"
#include "particledatabase.hpp"
#include "error.hpp"
#include "ibsimu.hpp"


using namespace std;


void test( int argc, char **argv )
{
    //Geometry geom( MODE_2D, Int3D(11,11,1), Vec3D(-0.05,-0.05,0.0), 0.01 );
    //Geometry geom( MODE_2D, Int3D(101,101,1), Vec3D(-0.05,-0.05,0.0), 0.001 );
    //Geometry geom( MODE_2D, Int3D(1001,1001,1), Vec3D(-0.05,-0.05,0.0), 0.0001 );
    Geometry geom( MODE_3D, Int3D(101,101,101), Vec3D(-0.05,-0.05,-0.05), 0.001 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 3, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 5, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 6, Bound(BOUND_DIRICHLET, 0.0) );
    geom.build_mesh();

    ScalarField epot( geom );
    MeshVectorField bfield;
    ScalarField scharge( geom );

    // Not needed
    //EpotProblem p;
    //p.construct( geom );
    //BiCGSTABSolver solver;
    //p.set_solver( solver );
    //p.solve( epot, scharge );

    EpotEfield efield( geom, epot );

    //ParticleDataBase2D pdb;
    ParticleDataBase3D pdb;
    pdb.set_thread_count( 1 );
    size_t n = 10;
    for( size_t a = 0; a < n; a++ ) {
	/*
	pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						      0.0, 5.0e4*(-0.5*n+a), 
						      -0.05, 1.0e5 ) );
	*/
	pdb.add_particle( 1.0, 1.0, 1.0, ParticleP3D( 0.0, 
						      0.0, 5.0e4*(-0.5*n+a), 
						      -0.05, 1.0e5,
						      -0.05, 5.0e4 ) );
    }
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    GeomPlotter geomplotter( geom );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.set_particle_div( 1 );
    geomplotter.plot_png( "particles_nofield.png" );

    /*
    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_scharge( &scharge );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
    */
}


