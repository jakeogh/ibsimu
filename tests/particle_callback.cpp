/*! \file particle_callback.cpp 
 *  \brief Test particle handler callback.
 *
 *  \test Test particle handler callback.
 */

#include <fstream>
#include <iomanip>
#include "gtkplotter.hpp"
#include "geometry.hpp"
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "func_solid.hpp"
#include "meshvectorfield.hpp"
#include "epot_efield.hpp"
#include "particles.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "ibsimutest.hpp"

#include "geomplotter.hpp"


using namespace std;

class THCallback : public TrajectoryHandlerCallback {
public:

    THCallback() {}

    virtual ~THCallback() {}

    virtual void operator()( ParticleBase *particle, ParticlePBase *xcur, ParticlePBase *xend ) const {
	ParticleP2D *pcur = (ParticleP2D *)( xcur );
	ParticleP2D *pend = (ParticleP2D *)( xend );
	//std::cout << "Callback: pcur = " << *pcur << "\n";
	if( (*pcur)[3] >= 0.05 ) {
	    //std::cout << "Collision\n";
	    *pend = *pcur;
	    particle->set_status( PARTICLE_COLL );
	}
    }
};


void test( int argc, char **argv )
{
    Geometry geom( MODE_2D, Int3D(71,71,1), Vec3D(0.0,0.0,0.0), 0.001 );

    geom.set_boundary( 1, Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,     0.0) );
    geom.set_boundary( 3, Bound(BOUND_DIRICHLET,   0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,  -1.0e3) );
    geom.build_mesh();

    EpotProblem p;
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    EpotEfield efield( geom, epot );
    MeshVectorField bfield;

    THCallback thc;

    ParticleDataBase2D pdb;
    pdb.set_trajectory_handler_callback( &thc );
    pdb.set_thread_count( 1 );

    pdb.add_particle( 0.0, 6.0, 19.9924, ParticleP2D( 0, 0.01, 3.0e4, 0.0, 3.0e4 ) );
    pdb.set_polyint( true );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    // Test for particle location
    Particle2D part = pdb.particle( 0 );
    if( fabs(part(0) - 4.24458e-7) > 1e-9 || fabs(part(1) - 0.0227322) > 1e-4 || 
	fabs(part(3) - 0.05) > 0.0001 || fabs(part(4) - 205591) > 100.0 )
	throw( ErrorTest( ERROR_LOCATION, "particle at incorrect location" ) );

    // Test last trajectory point
    size_t s = pdb.traj_size( 0 );
    ParticleP2D pp = pdb.trajectory_point( 0, s-1 );
    if( fabs(pp(0) - 4.24458e-7) > 1e-9 || fabs(pp(1) - 0.0227322) > 1e-4 || 
	fabs(pp(3) - 0.05) > 0.0001 || fabs(pp(4) - 205591) > 100.0 )
	throw( ErrorTest( ERROR_LOCATION, "particlepoint at incorrect location" ) );

    GeomPlotter gplotter( &geom );
    gplotter.set_size( 1024, 768 );
    gplotter.set_epot( &epot );
    gplotter.set_particle_div( 1 );
    gplotter.set_particle_database( &pdb );
    gplotter.plot_png( "particle_callback.png" );

    /*
    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &g );
    plotter.set_scharge( &scharge );
    plotter.set_epot( &epot );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
    */
}



