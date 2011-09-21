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
#include "random.hpp"
#include "gtkplotter.hpp"



using namespace std;


class SECallback : public TrajectoryEndCallback {

    Geometry &_geom;
    MTRandom  _rand;
    double    _k;
    double    _theta;
    
public:

    SECallback( Geometry &geom ) 
	: _geom(geom), _rand(2) {

	_k = 9.0;
	_theta = 0.5;
	_rand.set_transformation( 0, Cosine_Transformation() );
	_rand.set_transformation( 1, Gamma_Transformation( _k, _theta ) );
    }

    virtual ~SECallback() {

    }

    virtual void operator()( ParticleBase *particle, class ParticleDataBase *pdb ) const {

	Particle2D *p2d = (Particle2D *)( particle );
	Vec3D loc = p2d->location();
	//std::cout << "Particle end: " << p2d->x() << "\n";
	if( loc[0]*loc[0]+loc[1]*loc[1] <= 0.022*0.022 ) {
	    
	    // Get normal
	    Vec3D normal = _geom.surface_normal( loc );

	    // Adjust location off the surface
	    loc += 0.01*_geom.h()*normal;

	    // Launch 1000 particles
	    for( size_t a = 0; a < 1000; a++ ) {

		// Randomize velocity and direction
		double x[2];
		_rand.get( x );
		double angle = 0.5*M_PI*x[0];
		double mass = 1.0/1500.0;
		double speed = sqrt( 2.0*x[1]*CHARGE_E/(mass*MASS_U) );
		// Rotate normal according to cosine distribution
		Vec3D vel( cos(angle)*normal[0] - sin(angle)*normal[1],
			   sin(angle)*normal[0] + cos(angle)*normal[1],
			   0.0 );
		vel *= speed;

		ParticleDataBase2D *pdb2d = (ParticleDataBase2D *)( pdb );
		pdb2d->add_particle( 1.0, mass, 1.0, ParticleP2D( 0.0, 
								  loc[0], vel[0], 
								  loc[1], vel[1] ) );
	    }
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
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );

    SECallback secb( geom );
    pdb.set_trajectory_end_callback( &secb );

    // Crap def
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  0.1, 1.0e5,
						  0.1, 4.0e4 ) );

    // From left
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  -0.0499, 1.0e5,
						  0.0, 4.0e4 ) );
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  -0.0499, 1.0e5,
						  0.0, 1.0e4 ) );
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  -0.0499, 1.0e5,
						  0.0, -1.0e4 ) );
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  -0.0499, 1.0e5,
						  0.0, -4.0e4 ) );

    // From right
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  0.0499, -1.0e5,
						  0.0, 4.0e4 ) );
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  0.0499, -1.0e5,
						  0.0, 1.0e4 ) );
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  0.0499, -1.0e5,
						  0.0, -1.0e4 ) );
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  0.0499, -1.0e5,
						  0.0, -4.0e4 ) );

    // From top
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  0.0,  -2.0e4,
						  0.0499, -1.0e5 ) );
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  0.0,   2.0e4,
						  0.0499, -1.0e5 ) );

    // From bottom
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  0.0,  -2.0e4,
						  -0.0499, 1.0e5 ) );
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP2D( 0.0, 
						  0.0,   2.0e4,
						  -0.0499, 1.0e5 ) );
    
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    /*
    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_scharge( &scharge );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
    */
}
