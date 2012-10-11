/*! \file energy.cpp 
 *  \brief Test energy distribution with a beam in 2d system.
 *
 *  \test Test energy distribution with a beam in 2d system.
 */


#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "epot_gssolver.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_field.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "gtkplotter.hpp"


using namespace std;


void test( int argc, char **argv )
{    
    // 10x10 cm geometry with 0.25 cm mesh
    Geometry geom( MODE_2D, Int3D(41,41,1), Vec3D(0,0,0), 0.0025 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,  0.0) );
    geom.build_mesh();

    EpotGSSolver solver( geom );
    EpotField epot( geom );
    MeshScalarField scharge( geom );
    solver.solve( epot, scharge );

    EpotEfield efield( epot );
    MeshVectorField bfield;

    ParticleDataBase2D pdb;
    pdb.set_thread_count( 1 );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.add_2d_beam_with_energy( 1000, 50.0, 1.0, 1.0, 
				 3.0e3, 1.0, 0.0, 
				 0.0, 0.0, 
				 0.0, 0.009 );
    pdb.set_accuracy( 1.0e-9, 1.0e-9 );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    ParticleDiagPlotter diag( geom, pdb, AXIS_X, 0.0, PARTICLE_DIAG_PLOT_HISTO1D, 
			      DIAG_VX );
    diag.plot_png( "energy_velocity.png" );

    ParticleDiagPlotter diag2( geom, pdb, AXIS_X, 0.0, PARTICLE_DIAG_PLOT_HISTO1D, 
			      DIAG_EK );
    diag2.plot_png( "energy_energy.png" );

    if( false ) {
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.set_scharge( &scharge );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
    }
}
