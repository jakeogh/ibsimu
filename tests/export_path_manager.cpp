/*! \file export_path_manager.cpp 
 *  \brief Test path manager particle export
 *
 *  \test Test path manager particle export
 */


#include <iostream>
#include <cstdlib>
#include <iomanip>
#include "geometry.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "particles.hpp"
#include "particledatabase.hpp"
#include "error.hpp"
#include "ibsimu.hpp"

#include "gtkplotter.hpp"

using namespace std;


void test( int argc, char **argv )
{
    Geometry geom( MODE_3D, Int3D(51,51,52), Vec3D(-0.05,-0.05,0.0), 0.002 );
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

    EpotEfield efield( geom, epot );

    ParticleDataBase3D pdb;
    size_t n = 5;
    size_t m = 5;
    double vp = 2.0e6;
    double vt = 3.0e5;
    double vp_ref = 2.5e6;
    double Ep = 0.5*2.0*MASS_U*vp_ref*vp_ref/CHARGE_E;
    for( size_t a = 0; a < n; a++ ) {
	for( size_t b = 0; b < m; b++ ) {
	    pdb.add_particle( 1.0, 1.0, 2.0, ParticleP3D( 0.0, 
							  0.0, vt*(-0.5*(n-1)+a) / (0.5*(n-1)), 
							  0.0, vt*(-0.5*(m-1)+b) / (0.5*(m-1)), 
							  0.0, vp ) );
	}
    }
    pdb.iterate_trajectories( scharge, efield, bfield, geom );
    pdb.export_path_manager_data( "export_path_manager.txt", Ep, 1.0, 2.0, 
				  Vec3D(0,0,0.1),
				  Vec3D(1,0,0),
				  Vec3D(0,1,0) );
    //pdb.debug_print( cout );

    /*

      The test file contains 25 particles with charge of 1 electron
      charge and mass of 2 atomic mass units. Each of the particles is
      carrying 1 A current.

      The particles originate from origin (0,0,0) and have velocity of
      2e6 m/s in z-direction and a flat spread of velocity in
      transverse direction (x and y) from -3e5 m/s to +3e5 m/s. The x
      and y velocities are therefore -3.0e5, -1.5e5, 0.0, +1.5e5 and
      +3.0e5 in units m/s. The output export is made on z=0.1 m plane.

      The x' and y' angles should range from 3.0e5/2e6 = 0.15 rad to
      -3.0e5/2e6 = -0.15 rad. The x and y should range from 0.15*0.1m
      = 0.015 m to -0.15*0.1m = -0.015 m.

      The reference momentum is set for velocity 2.5e6 m/s. All the
      particles have p_z velocity 2.0e6 so th dp/p_ref should be
      -0.2 for all particles. 

      Phase is randomized between -pi and pi.

     */

    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_scharge( &scharge );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
}


