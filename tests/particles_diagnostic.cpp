/*! \file particles_diagnostic.cpp 
 *  \test Test particle diagnostic tool.
 *
 */


#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "epot_efield.hpp"
#include "particles.hpp"
#include "particledatabase.hpp"
#include "error.hpp"
#include "ibsimu.hpp"


using namespace std;


#define ERROR()							    \
    {								    \
	cout << "Error at " << __FILE__ << ":" << __LINE__ << "\n"; \
	exit( 1 );						    \
    }


void test( void )
{
    //verbose_output = 1;

    Geometry geom( MODE_2D, Int3D(11,11,1), Vec3D(0,0,0), 0.01 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET, 1000.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,    0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,      0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,      0.0) );
    geom.build_mesh();

    EpotProblem p;
    p.construct( geom );

    ScalarField epot( geom );
    VectorField bfield;
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    EpotEfield efield( geom, epot );

    ParticleDataBase2D pdb;
    pdb.set_thread_count( 1 );
    pdb.add_particle( 0.0, 1.0, 1.0, ParticleP2D( 0, 0.01, 0.0, 0.0, 1e5 ) );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );
    //pdb.debug_print();

    std::vector<trajectory_diagnostic_e> diagnostics;
    diagnostics.push_back( DIAG_T  );
    diagnostics.push_back( DIAG_X  );
    diagnostics.push_back( DIAG_Y  );
    diagnostics.push_back( DIAG_VX  );
    diagnostics.push_back( DIAG_VY  );
    diagnostics.push_back( DIAG_XP );
    diagnostics.push_back( DIAG_YP );
    TrajectoryDiagnosticData tdata;
    pdb.trajectories_at_plane( tdata, AXIS_X, 0.08, diagnostics );

    // Check results

    // t
    double t = sqrt(((0.08-0.01)*2.0*MASS_U)/(CHARGE_E*10000.0));
    if( fabs(tdata(0,0)-t)/t > 0.001 )
	ERROR();

    // x
    if( fabs(tdata(0,1)-0.08)/0.08 > 0.001 )
	ERROR();

    // y
    if( fabs(tdata(0,2)-1e5*t)/(1e5*t) > 0.001 )
	ERROR();

    // vx
    double vx = CHARGE_E*10000.0/MASS_U*t;
    if( fabs(tdata(0,3)-vx)/vx > 0.001 )
	ERROR();

    // vy
    double vy = 1e5;
    if( fabs(tdata(0,4)-vy)/vy > 0.001 )
	ERROR();

    // xp
    if( fabs(tdata(0,5)-(vx/vx))/(vx/vx) > 0.001 )
	ERROR();

    // yp
    if( fabs(tdata(0,6)-(vy/vx))/(vy/vx) > 0.001 )
	ERROR();

    GeomPlotter geomplotter( &geom );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.plot_png( "particles_diagnostic.png" );
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



