/*! \file nsimp_plasmacyl.cpp
 *  \brief Test with a simple negative ion plasma in cylindrical 
 *  electrode configuration.
 *
 *  \test Test with a simple negative ion plasma in cylindrical 
 *  electrode configuration.
 */


#include <cstdlib>
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
#include "ibsimu.hpp"
#include "error.hpp"
#include "particlediagplotter.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"


using namespace std;


bool solid1( double x, double y, double z )
{
    return( x <= 0.004 && y >= 0.0030 && 
	    y >= 4.8*x - 0.003 );
}


bool solid2( double x, double y, double z )
{
    return( (x >= 0.025 || y <= 3.25*x - 0.06925) &&
	    x >= 0.0231 && y >= 0.0055 );
}


bool initial_plasma( double x, double y, double z )
{
    return( x <= 0.00055 );
}


void test( int argc, char **argv )
{
    // 30x20 mm geometry with 0.05 mm mesh size
    //Geometry geom( MODE_CYL, Int3D(601,401,1), Vec3D(0,0,0), 0.00005 );
    Geometry geom( MODE_CYL, Int3D(301,201,1), Vec3D(0,0,0), 0.0001 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0 ) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, +4.0e3) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0)  );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, +4.0e3) );
    geom.build_mesh();

    EpotProblem p;
    p.set_nsimp_initial_plasma( initial_plasma );
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );

    VectorField bfield;
    EpotEfield efield( geom, epot );
    field_extrpl_e efldextrpl[6] = {FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				    FIELD_MIRROR,      FIELD_EXTRAPOLATE,
				    FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBaseCyl pdb;
    pdb.set_thread_count( 4 );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );

    for( size_t i = 0; i < 2; i++ ) {

	if( i == 1 ) {
	    std::vector<double> Ei;
	    std::vector<double> rhoi;
	    double rhop = pdb.get_rhosum();
	    p.set_nsimp_plasma( rhop, 2.5, rhoi, Ei );
	    p.construct( geom );
	    solver.set_newton_imax( 4 );
	}

	p.solve( epot, scharge );

	pdb.clear();
	// H-
	// 1 mA total -> 35.37 A/m2
	pdb.add_2d_beam_with_energy( 5000, -35.37, -1.0, 1.0, 
				     5.0, 0.0, 0.5, 
				     0.0, 0.0, 
				     0.0, 0.003 );
	// e-
	// 20 mA total -> 707.4 A/m2
	pdb.add_2d_beam_with_energy( 5000, -707.4, -1.0, 1.0/1836.15, 
				     5.0, 0.0, 0.5, 
				     0.0, 0.0, 
				     0.0, 0.003 );
	pdb.iterate_trajectories( scharge, efield, bfield, geom );

	/*
	if( i == 4 ) {
	    GTKPlotter plotter( argc, argv );
	    plotter.set_geometry( &geom );
	    plotter.set_epot( &epot );
	    plotter.set_scharge( &scharge );
	    plotter.set_particledatabase( &pdb );
	    plotter.new_geometry_plot_window();
	    plotter.run();
	}
	*/

	GeomPlotter gplotter( &geom );
	gplotter.set_size( 1024, 768 );
	gplotter.set_epot( &epot );
	std::vector<double> eqlines;
	eqlines.push_back( -8.0 );
	eqlines.push_back( -4.0 );
	eqlines.push_back( 0.0 );
	eqlines.push_back( +4.0 );
	gplotter.set_eqlines_manual( eqlines );
	gplotter.set_particle_database( &pdb );
	string filename = (string)"nsimp_plasmacyl_" + to_string(i) + ".png";
	gplotter.plot_png( filename );

    }
}

