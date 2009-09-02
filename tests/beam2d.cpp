/*! \file beam2d.cpp 
 *  \brief Test with a beam in 2d system.
 *
 *  \test Test with a beam in 2d system. Space charge on axis is checked.
 *
 */


#include <fstream>
#include <iomanip>

#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "vectorfield.hpp"
#include "particledatabase.hpp"
#include "verbose.hpp"
#include "error.hpp"
#include "gtkplotter.hpp"


using namespace std;


void test( int *argc, char ***argv )
{
    verbose_output = 1;
    
    // 10x10 cm geometry with 0.25 cm mesh
    Geometry g( MODE_2D, Int3D(41,41,1), Vec3D(0,0,0), 0.0025 );
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

    ParticleDataBase2D pdb;
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
    ofstream ostr( "beam2d_map.dat" );
    ofstream ostr2( "beam2d.dat" );
    ostr << "# "
         << setw(12) << "x (m)" << " " 
         << setw(14) << "y (m)" << " " 
         << setw(14) << "potential (V)" << " "
         << setw(14) << "scharge (C/m3)" << "\n";
    ostr2 << "# "
         << setw(12) << "y (m)" << " " 
         << setw(14) << "potential (V)" << " "
         << setw(14) << "scharge (C/m3)" << "\n";
    bool err = false;
    for( int j = 0; j < g.size(1); j++ ) {
	Vec3D x( 0.05, g.h()*j, 0.0  );
	ostr2 << setw(14) << x[1] << " "
	      << setw(14) << epot( x ) << " "
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

    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &g );
    plotter.set_epot( &epot );
    plotter.set_scharge( &scharge );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();

    if( err ) {
	std::cout << "Error: calculated space charge differs from theory\n";
	exit( 1 );
    }
}


int main( int argc, char **argv )
{
    try {
	test( &argc, &argv );
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

