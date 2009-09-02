/*! \file particles_bfieldcyl.cpp 
 *  \brief Test particle iterator in constant magnetic field in 
 *  cylindrical geometry.
 *
 *  \test  Test particle iterator in constant magnetic field 
 *  and no electric field. %Particle trajectory is checked. The
 *  particle trajectory should be circular with radius
 *  \f[ r = \frac{mv}{qB} = \frac{1\mathrm{~u} \cdot 2\cdot10^6}
    {1\mathrm{~e} \cdot 1\mathrm{~T}} = 0.020728544\mathrm{~m} \f]
 *
 */


#include <fstream>
#include <iomanip>
#include "fileplot.hpp"
#include "geometry.hpp"
#include "bicgstab_solver.hpp"
#include "vectorfield.hpp"
#include "epot_problem.hpp"
#include "epot_efield.hpp"
#include "particles.hpp"
#include "error.hpp"
#include "verbose.hpp"


using namespace std;


void test( void )
{
    verbose_output = 1;

    Geometry g( MODE_CYL, Int3D(21,11,1), Vec3D(-0.05,0.0,0.0), 0.01 );
    g.set_boundary( 1, Bound(BOUND_DIRICHLET,    0.0) );
    g.set_boundary( 2, Bound(BOUND_DIRICHLET,    0.0) );
    g.set_boundary( 3, Bound(BOUND_DIRICHLET,    0.0) );
    g.set_boundary( 4, Bound(BOUND_DIRICHLET,    0.0) );
    g.build_mesh();
    //g.debug_print();

    EpotProblem p;
    p.construct( g );

    ScalarField epot( g );
    ScalarField scharge( g );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    EpotEfield efield( g, epot );

    // Make magnetic field
    bool fout[3] = {true,false,false};
    VectorField bfield( MODE_CYL, fout, Int3D(2,2,1), Vec3D(-0.05,0.0,0.0), 0.2 );
    bfield.set( 0, 0, Vec3D(1,0,0) );
    bfield.set( 1, 0, Vec3D(1,0,0) );
    bfield.set( 0, 1, Vec3D(1,0,0) );
    bfield.set( 1, 1, Vec3D(1,0,0) );

    ParticleDataBaseCyl pdb;
    pdb.set_thread_count( 1 );
    pdb.set_max_steps( 1000 );
    //pdb.set_max_time( 260e-9 );
    pdb.add_particle( 0.0, 1.0, 1.0, ParticlePCyl( 0, 
						   -0.045, 2e6, 
						   0.020728544449, 0.0, 
						   -2e6/0.020728544449 ) );
    pdb.add_particle( 0.0, 1.0, 1.0, ParticlePCyl( 0, 
						   -0.045, 2e6, 
						   0.04, 0.0, 
						   -2e6/0.04 ) );
    pdb.add_particle( 0.0, 1.0, 1.0, ParticlePCyl( 0, 
						   -0.045, 2e6, 
						   0.01, 0.0, 
						   -2e6/0.01 ) );
    pdb.iterate_trajectories( scharge, efield, bfield, g );
    //pdb.debug_print();

    // Check particle trajectory points
    bool err = false;
    ofstream ostr( "particles_bfieldcyl.dat" );
    ostr << "# "
	 << setw(12) << "Time (s)" << " "
	 << setw(14) << "x (m)" << " "
	 << setw(14) << "r (m)" << "\n";
    ParticleCyl &prt = pdb.particle(0);
    for( uint32_t b = 0; b < prt.traj_size(); b++ ) {
	double t = prt.traj(b)(0);
	double x = prt.traj(b)(1);
	double r = prt.traj(b)(3);
	ostr << setw(14) << t << " "
	     << setw(14) << x << " "
	     << setw(14) << r << "\n";
	if( fabs( r - 0.020728544449 ) > 1e-5 ) {
	    std::cout << "Error: trajectory differs from theory\n";
	    exit( 1 );
	}
    }

    GeomPlotter geomplotter;
    geomplotter.set_geometry( g );
    geomplotter.set_epot( epot );
    geomplotter.set_particledatabase( pdb );
    pngplot( &geomplotter, "particles_bfieldcyl.png" );

    ostr.close();
    if( err ) {
	std::cout << "Error: trajectory differs from theory\n";
	exit( 1 );
    }
}


int main( void )
{
    try {
	test();
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

