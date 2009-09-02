/*! \file plasmacyl.cpp 
 *  \brief Test with a plasma in cylindrical electrode configuration.
 *
 *  \test Test with a plasma in cylindrical electrode configuration.
 *
 */


#include <sstream>
#include <fstream>
#include <iomanip>
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "particledatabase.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "vectorfield.hpp"
#include "verbose.hpp"
#include "error.hpp"

#include "gtkplotter.hpp"
#include "frame.hpp"
#include "solidplot.hpp"
#include "particleplot.hpp"
#include "eqpotplot.hpp"

using namespace std;


bool solid1( double x, double y, double z )
{
    return( x <= 0.00187 && y >= 0.00054 && y >= 2.28*x - 0.0010 &&
	    (x >= 0.00054 || y >= 0.0015) );
}


bool solid2( double x, double y, double z )
{
    return( x >= 0.0095 && y >= 0.0023333 && y >= 0.01283 - x );
}


void test( int *argc, char ***argv )
{
    verbose_output = 1;
    
    // 12x7 mm geometry with 0.05 mm mesh size
    Geometry g( MODE_CYL, Int3D(241,141,1), Vec3D(0,0,0), 0.00005 );
    Solid *s1 = new FuncSolid( solid1 );
    g.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    g.set_solid( 8, s2 );
    g.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0 ) );
    g.set_boundary( 2, Bound(BOUND_DIRICHLET, -8.0e3) );
    g.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0)  );
    g.set_boundary( 8, Bound(BOUND_DIRICHLET, -8.0e3) );
    g.build_mesh();

    EpotProblem p;
    p.set_initial_plasma( 5.0, 0.00055 );
    p.construct( g );

    ScalarField epot( g );
    ScalarField scharge( g );

    BiCGSTABSolver solver;
    p.set_solver( solver );

    VectorField bfield;
    EpotEfield efield( g, epot );
    efield_extrpl_e efldextrpl[6] = {EFIELD_EXTRAPOLATE, EFIELD_EXTRAPOLATE, 
				     EFIELD_MIRROR,EFIELD_EXTRAPOLATE,
				     EFIELD_EXTRAPOLATE, EFIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBaseCyl pdb;
    pdb.set_thread_count( 4 );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );

    /*
    GeomPlotter geomplotter;
    geomplotter.set_size( 1024, 768 );
    //geomplotter.set_range( 0.09, 0.12000001, 0.02, 0.04 );
    //geomplotter.set_meshlines( true );
    geomplotter.set_geometry( g );
    geomplotter.set_epot( epot );
    std::vector<double> eqlines;
    eqlines.push_back( -8.0 );
    eqlines.push_back( -4.0 );
    eqlines.push_back( 0.0 );
    eqlines.push_back( +4.0 );
    geomplotter.set_manual_eqlines( eqlines );
    geomplotter.set_particledatabase( pdb );
    */

    for( size_t i = 0; i < 1; i++ ) {

	if( i == 1 ) {
	    double rhoe = pdb.get_rhosum();
	    p.set_pexp_plasma( -rhoe, 5.0, 5.0 );
	    p.construct( g );
	}

	p.solve( epot, scharge );
	pdb.clear();
	pdb.add_2d_beam_with_energy( 5000, 600.0, 1.0, 1.0, 
				     5.0, 0.0, 0.5, 
				     0.0, 0.0, 
				     0.0, 0.0015 );
	pdb.iterate_trajectories( scharge, efield, bfield, g );
	//stringstream ss;
	//ss << "plasmacyl_" << i << ".png";
	//pngplot( &geomplotter, ss.str() );
    }

    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &g );
    plotter.set_epot( &epot );
    plotter.set_scharge( &scharge );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
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

