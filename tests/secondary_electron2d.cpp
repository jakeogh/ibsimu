/*! \file secondary_electron2d.cpp 
 *  \brief Test secondary electron emission in 2d.
 *
 *  \test  Test secondary electron emission in 2d.
 */


#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "particledatabase.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "gtkplotter.hpp"


using namespace std;


class SECallback : public TrajectoryEndCallback {
public:

    virtual ~SECallback() {}

    virtual void operator()( ParticleBase *particle, class ParticleDataBase *pdb ) const {

	Particle2D *p2d = (Particle2D *)( particle );
	Vec3D loc = p2d->location();
	std::cout << "Particle end: " << p2d->x() << "\n";
	if( loc[0]*loc[0]+loc[1]*loc[1] <= 0.022*0.022 ) {
	
	    ParticleDataBase2D *pdb2d = (ParticleDataBase2D *)( pdb );
	    pdb2d->add_particle( 1.0, 1.0, 1.0, ParticleP2D( (*p2d)[0], 
							     (*p2d)[1], 0.0, 
							     (*p2d)[3], 1.0e6 ) );
	}
    }
};


bool solid1( double x, double y, double z )
{
    return( x*x+y*y <= 0.02*0.02 );
}


void test( int argc, char **argv )
{    
    Geometry geom( MODE_2D, Int3D(101,101,1), Vec3D(-0.05,-0.05,0), 0.001 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 3, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
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

    ParticleDataBase2D pdb;
    pdb.set_thread_count( 1 );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );

    SECallback secb;
    pdb.set_trajectory_end_callback( &secb );

    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  -0.05, 1.0e5,
						  0.0, 1.0e4 ) );
    
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_scharge( &scharge );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
}
