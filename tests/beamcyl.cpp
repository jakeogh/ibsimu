/*! \file beamcyl.cpp 
 *  \brief Test with a beam in cylindrical system.
 *
 *  \test Test with a beam in cylindrical system. Space charge on axis is checked.
 *
 */


#include <fstream>
#include <iomanip>
#include "fileplot.hpp"
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "verbose.hpp"
#include "error.hpp"


using namespace std;


void test( void )
{
    verbose_output = 1;
    
    // 10x10 cm geometry with 0.25 cm mesh
    Geometry g( MODE_CYL, Int3D(41,41,1), Vec3D(0,0,0), 0.0025 );
    g.set_boundary( 1, Bound(BOUND_NEUMANN,  0.0) );
    g.set_boundary( 2, Bound(BOUND_NEUMANN,  0.0) );
    g.set_boundary( 3, Bound(BOUND_NEUMANN,  0.0) );
    g.set_boundary( 4, Bound(BOUND_DIRICHLET,  0.0) );
    g.build_mesh();

    EpotProblem p;
    p.construct( g );

    ScalarField epot( g );
    ScalarField scharge( g );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    EpotEfield efield( g, epot );
    VectorField bfield;

    ParticleDataBaseCyl pdb;
    pdb.set_thread_count( 1 );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.add_2d_beam_with_energy( 1000, 50.0, 1.0, 1.0, 
				 3.0e3, 0.0, 0.0, 
				 0.0, 0.0, 
				 0.0, 0.009 );
    pdb.iterate_trajectories( scharge, efield, bfield, g );
    p.solve( epot, scharge );

    // Write calculated space charge and epot to file
    ofstream ostr( "beamcyl_map.dat" );
    ofstream ostr2( "beamcyl.dat" );
    ostr << "# "
         << setw(12) << "x (m)" << " " 
         << setw(14) << "r (m)" << " " 
         << setw(14) << "potential (V)" << " "
         << setw(14) << "scharge (C/m3)" << "\n";
    ostr2 << "# "
         << setw(12) << "r (m)" << " " 
         << setw(14) << "potential (V)" << " "
         << setw(12) << "E_x (V/m)" << " "
         << setw(12) << "E_r (V/m)" << " "
         << setw(12) << "E_z (V/m)" << " "
         << setw(14) << "scharge (C/m3)" << "\n";
    bool err = false;
    for( int j = 0; j < g.size(1); j++ ) {
	Vec3D x( 0.05, g.h()*j, 0.0  );
	ostr2 << setw(14) << x[1] << " "
	      << setw(14) << epot( x ) << " "
	      << setw(14) << efield( x ) << " "
	      << setw(14) << scharge( x ) << "\n";
	for( int i = 0; i < g.size(0); i++ ) {
	    Vec3D x( g.h()*i, g.h()*j, 0.0  );
	    if( j == 0 && fabs(scharge(x)-6.57148899428e-5) > 1e-9 )
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

    GeomPlotter geomplotter;
    geomplotter.set_size( 640, 480 );
    //geomplotter.set_range( 0.09, 0.12000001, 0.02, 0.04 );
    geomplotter.set_meshlines( true );
    geomplotter.set_geometry( g );
    geomplotter.set_epot( epot );
    geomplotter.set_particledatabase( pdb );
    pngplot( &geomplotter, "beamcyl.png" );

    if( err ) {
	std::cout << "Error: calculated space charge differs from theory\n";
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

