/*! \file particles_mirror.cpp
 *  \brief Test particle mirroring
 *
 *  \test Test particle mirroring
 *
 *  The test case has a 1000 V potential difference between
 *  x=-0.05 m and x=0.05 m planes. The electric field is therefore
 *  -1000 V / 0.1 m = -10000 V/m. The particle trajectories should
 *  follow trajectories defined by constant acceleration motion
 *  \f[ x = x_0 + v_0 t + \frac{1}{2} a_x t^2, \f]
 *  where the acceleration is
 *  \f[ a_x = \frac{F_x}{m} = \frac{Eq}{m}. \f]
 *  In this test case the particles have q = 1e, m = 1u. Therefore 
 *  \f$a_x = -9.64853082148e11\mathrm{~m/s}^2 \f$
 */


#include <iostream>
#include <iomanip>
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "epot_efield.hpp"
#include "particles.hpp"
#include "error.hpp"
#include "ibsimu.hpp"


using namespace std;


const double ax = -9.64853082148e11;
double maxerr = 0.0;
bool err = false;


void check_particle( Particle2D &p, double x0, double vx, double y0, double vy )
{
    for( uint32_t b = 0; b < p.traj_size(); b++ ) {
	double t = p.traj(b)(0);
	//std::cout << t << "\n";
	if( t < 3.21935897726e-7 - 0.1e-7 ) {
	    double x = x0 + vx * t + 0.5*ax*t*t;
	    double y = y0 + vy * t;
	    //std::cout << x << "\t" << y << "\t" << p.traj(b)(1) << "\t" << p.traj(b)(3) << "\n";
	    double xerr = fabs( p.traj(b)(1) - x );
	    double yerr = fabs( p.traj(b)(3) - y );
	    if( xerr > 5e-5 || yerr > 5e-5 )
		err = true;
	    //std::cout << "xerr = " << xerr << "\n";
	    //std::cout << "yerr = " << yerr << "\n";
	} else if( t > 3.21935897726e-7 + 0.1e-7 ) {
	    double t2 = 3.21935897726e-7;
	    double new_x0 = x0 + vx * t2 + 0.5*ax*t2*t2;
	    double new_y0 = y0 + vy * t2;
	    double new_vx = -vx - ax*t2;
	    t = t - t2;
	    double x = new_x0 + new_vx * t + 0.5*ax*t*t;
	    double y = new_y0 + vy * t;
	    double xerr = fabs( p.traj(b)(1) - x );
	    double yerr = fabs( p.traj(b)(3) - y );
	    //std::cout << x << "\t" << y << "\t" << p.traj(b)(1) << "\t" << p.traj(b)(3) << "\n";
	    if( xerr > 5e-5 || yerr > 5e-5 ) 
		err = true;
	    //std::cout << "xerr = " << xerr << "\n";
	    //std::cout << "yerr = " << yerr << "\n";
	}	    
    }
}


void test( int argc, char **argv )
{
    Geometry geom( MODE_2D, Int3D(11,11,1), Vec3D(-0.05,-0.05,0.0), 0.01 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 1000.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,      0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,      0.0) );
    geom.build_mesh();
    //g.debug_print();

    EpotProblem p;
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );
    VectorField bfield;

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    EpotEfield efield( geom, epot );
    field_extrpl_e efldextrpl[6] = { FIELD_MIRROR, FIELD_MIRROR, 
				     FIELD_MIRROR, FIELD_MIRROR, 
				     FIELD_MIRROR, FIELD_MIRROR };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase2D pdb;
    pdb.set_thread_count( 1 );
    bool pmirror[6] = { true, false, false, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.add_particle( 0.0, 1.0, 1.0, ParticleP2D( 0.0, 0.0, 0.0, -0.04, 1e5 ) );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );
    //pdb.debug_print();

    // Check particle trajectory points
    check_particle( pdb.particle(0), 0.0, 0.0, -0.04, 1e5 );

    GeomPlotter geomplotter( &geom );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.plot_png( "particles_mirror.png" );

    if( err ) {
	std::cout << "Error: trajectory differs from theory before mirroring\n";
	exit( 1 );
    }
}


