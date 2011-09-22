/*! \file secondary_electron3d.cpp 
 *  \brief Test secondary electron emission in 3d.
 *
 *  \test  Test secondary electron emission in 3d.
 */


#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
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
#include "geomplotter.hpp"
#include "random.hpp"


using namespace std;


double maxerr = 0.0;


class SECallback : public TrajectoryEndCallback {

    Geometry &_geom;
    MTRandom  _rand;
    double    _k;
    double    _theta;
    
public:

    SECallback( Geometry &geom ) 
	: _geom(geom), _rand(3) {

	_k = 9.0;
	_theta = 0.5;
	_rand.set_transformation( 0, Gamma_Transformation( _k, _theta ) );
	//_rand.set_transformation( 1, Cosine_Transformation() );
    }

    virtual ~SECallback() {}

    virtual void operator()( ParticleBase *particle, class ParticleDataBase *pdb ) const {

	Particle3D *p3d = (Particle3D *)( particle );
	Vec3D loc = p3d->location();
	//std::cout << "Particle end: " << p3d->x() << "\n";
	if( loc[0]*loc[0] + loc[1]*loc[1] + loc[2]*loc[2] <= 0.022*0.022 ) {

	    // Get normal
	    Vec3D normal = _geom.surface_normal( loc );

	    // Adjust location off the surface
	    loc += 0.01*_geom.h()*normal;

	    // Launch 50000 particles
	    for( size_t a = 0; a < 50000; a++ ) {

		// Randomize velocity and direction
		double x[2];
		_rand.get( x );
		double mass = 1.0/1500.0;
		double speed = sqrt( 2.0*x[0]*CHARGE_E/(mass*MASS_U) );

		// Find tangents
		Vec3D tang1 = normal.arb_perpendicular();
		Vec3D tang2 = cross( normal, tang1 );
		tang1.normalize();
		tang2.normalize();
		//std::cout << "normal = " << normal << "\n";
		//std::cout << "tang1 = " << tang1 << "\n";
		//std::cout << "tang2 = " << tang2 << "\n";

		// Build velp in natural coordinates
		double azm_angle = 2.0*M_PI*x[1];
		double pol_angle = asin( sqrt(x[2]) );
		Vec3D velp( speed*cos(pol_angle),  speed*sin(pol_angle), 0.0 );
		Transformation t;
		t.rotate_x( azm_angle );
		velp = t.transform_vector( velp );

		// Convert to surface coordinates
		Vec3D vel = velp[0]*normal + velp[1]*tang1 + velp[2]*tang2;

		ParticleDataBase3D *pdb3d = (ParticleDataBase3D *)( pdb );
		pdb3d->add_particle( 1.0, mass, 1.0, ParticleP3D( 0.0, 
								  loc[0], vel[0], 
								  loc[1], vel[1], 
								  loc[2], vel[2] ) );
	    }
	}
    }
};


bool solid1( double x, double y, double z )
{
    return( x*x+y*y+z*z <= 0.02*0.02 );
}


void test( int argc, char **argv )
{    
    Geometry geom( MODE_3D, Int3D(101,101,101), Vec3D(-0.05,-0.05,-0.05), 0.001 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 3, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 5, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 6, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    geom.build_mesh();

    ibsimu.set_thread_count( 1 );

    EpotProblem p;
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    EpotEfield efield( geom, epot );
    MeshVectorField bfield;

    ParticleDataBase3D pdb;
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );

    SECallback secb( geom );
    pdb.set_trajectory_end_callback( &secb );

    // From left
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP3D( 0.0, 
						  -0.0499, 1.0e5,
						  0.0, 3.0e4,
						  0.0, 2.0e4 ) );

    // From right
    pdb.add_particle( 1.0, 1.0, 1.0, ParticleP3D( 0.0, 
						  0.0499, -1.0e5,
						  0.0, 1.0e4,
						  0.0, -4.0e4 ) );


    
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    //std::cout << "Maxerr = " << maxerr << "\n";

    if( false ) {
	ScalarField tdens( geom );
	pdb.build_trajectory_density_field( tdens );
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.set_scharge( &scharge );
	plotter.set_trajdens( &tdens );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
    }

    if( false ) {
	ScalarField tdens( geom );
	pdb.build_trajectory_density_field( tdens );
	GeomPlotter gplotter( &geom );
	gplotter.set_trajdens( &tdens );
	gplotter.set_fieldgraph_plot( FIELD_TRAJDENS );
	gplotter.set_fieldgraph_logscale( true );
	gplotter.set_size( 800, 600 );
	gplotter.set_font_size( 16 );
	
	for( int level = 0; level < geom.size(2); level++ ) {
	    gplotter.set_view( VIEW_XY, level );
	    stringstream ss_filename;
	    ss_filename <<  "plot_" << setfill('0') << setw(3) << level << ".png";
	    gplotter.plot_png( ss_filename.str() );
	}
    }
}
