/*! \file beamcyl.cpp 
 *  \brief Test with a beam in cylindrical system. Space 
 *  charge on axis is checked.
 *
 *  \test Test with a beam in cylindrical system. Space 
 *  charge on axis is checked.
 */


#include <fstream>
#include <iomanip>
#include "geomplotter.hpp"
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "particlediagplotter.hpp"



using namespace std;


void test( int argc, char **argv )
{
    // 10x10 cm geometry with 0.25 cm mesh
    Geometry geom( MODE_CYL, Int3D(41,41,1), Vec3D(0,0,0), 0.0025 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,  0.0) );
    geom.build_mesh();

    EpotProblem p;
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    EpotEfield efield( geom, epot );
    VectorField bfield;

    ParticleDataBaseCyl pdb;
    pdb.set_thread_count( 1 );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    // Neon (m=20) 6+ beam of J = 45 A/m2 with r = 9 mm and R = 3 keV.
    // Space charge density is rho = J / v = J * sqrt(m/(2E)) = 2.64497329883e-4 C/m^3
    pdb.add_2d_beam_with_energy( 1000, 45.0, 6.0, 20.0, 
				 3.0e3, 0.0, 0.0, 
				 0.0, 0.0, 
				 0.0, 0.009 );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );
    p.solve( epot, scharge );

    // Write calculated space charge and epot to file
    ofstream ostr( "beamcyl_map.dat" );
    ostr << "# "
         << setw(12) << "x (m)" << " " 
         << setw(14) << "r (m)" << " " 
         << setw(14) << "potential (V)" << " "
         << setw(14) << "scharge (C/m3)" << "\n";
    bool err = false;
    for( int j = 0; j < geom.size(1); j++ ) {
	for( int i = 0; i < geom.size(0); i++ ) {
	    Vec3D x( geom.h()*i, geom.h()*j, 0.0  );
	    if( j == 0 && fabs(scharge(x)-2.64497329883e-4) > 1e-9 )
		err = true;
	    ostr << setw(14) << x[0] << " "
		 << setw(14) << x[1] << " "
		 << setw(14) << epot( x ) << " "
		 << setw(14) << scharge( x ) << "\n";
	}
	ostr << "\n";
    }
    ostr << "\n\n";
    ostr.close();

    ostr.open( "beamcyl.dat" );
    ostr << "# "
         << setw(12) << "r (m)" << " " 
         << setw(14) << "potential (V)" << " "
         << setw(12) << "E_x (V/m)" << " "
         << setw(12) << "E_r (V/m)" << " "
         << setw(12) << "E_z (V/m)" << " "
         << setw(14) << "scharge (C/m3)" << "\n";
    for( double y = 0.0; y < 0.1; y += 0.0005 ) {
	Vec3D x( 0.05, y, 0.0 );
	ostr << setw(14) << y << " "
	     << setw(14) << epot( x ) << " "
	     << setw(14) << efield( x ) << " "
	     << setw(14) << scharge( x ) << "\n";
    }
    ostr << "\n\n";
    ostr.close();

    // Do current density plot
    ParticleDiagPlotter pplotter( &geom, &pdb, AXIS_X, 0.05, PARTICLE_DIAG_PLOT_HISTO1D, DIAG_R );
    pplotter.set_font_size( 18 );
    pplotter.set_histogram_n( 101 );
    pplotter.set_size( 640, 640 );
    pplotter.export_data( "beamcyl_curr.dat" );
    pplotter.plot_png( "beamcyl_curr.png" );

    GeomPlotter geomplotter( &geom );
    geomplotter.set_size( 640, 480 );
    geomplotter.set_mesh( true );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.plot_png( "beamcyl.png" );

    if( err ) {
	std::cout << "Error: calculated space charge differs from theory\n";
	exit( 1 );
    }
}
