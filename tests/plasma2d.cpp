/*! \file plasma2d.cpp 
 *  \brief Test with a plasma in 2d electrode configuration.
 *
 *  \test Test with a plasma in 2d electrode configuration.
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
#include "geomplotter.hpp"
#include "particlediagplotter.hpp"
#include "fielddiagplotter.hpp"


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
    Geometry geom( MODE_2D, Int3D(241,141,1), Vec3D(0,0,0), 0.00005 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0 ) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, -8.0e3) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0)  );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, -8.0e3) );
    geom.build_mesh();

    EpotProblem p;
    p.set_initial_plasma( 5.0, 0.00055 );
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );

    VectorField bfield;
    EpotEfield efield( geom, epot );
    efield_extrpl_e efldextrpl[6] = {EFIELD_EXTRAPOLATE, EFIELD_EXTRAPOLATE, 
				     EFIELD_MIRROR,EFIELD_EXTRAPOLATE,
				     EFIELD_EXTRAPOLATE, EFIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase2D pdb;
    pdb.set_thread_count( 4 );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );

    for( size_t i = 0; i < 2; i++ ) {

	if( i == 1 ) {
	    double rhoe = pdb.get_rhosum();
	    p.set_pexp_plasma( -rhoe, 5.0, 5.0 );
	    p.construct( geom );
	}

	p.solve( epot, scharge );
	pdb.clear();
	pdb.add_2d_beam_with_energy( 30000, 600.0, 1.0, 1.0, 
				     5.0, 0.0, 0.5, 
				     0.0, 0.0, 
				     0.0, 0.0015 );
	pdb.iterate_trajectories( scharge, efield, bfield, geom );
    }

    /*
    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_scharge( &scharge );
    plotter.set_epot( &epot );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
    */

    GeomPlotter gplotter( &geom );
    gplotter.set_scharge( &scharge );
    gplotter.set_scharge_field( true );
    gplotter.set_epot( &epot );
    gplotter.set_particle_database( &pdb );
    gplotter.set_particle_div( 0 );
    std::vector<double> pot;
    pot.push_back( +2 );
    pot.push_back( 0 );
    pot.push_back( -2 );
    pot.push_back( -20 );
    gplotter.set_eqlines_manual( pot );
    gplotter.set_font_size( 15 );
    gplotter.plot_png( "plasma2d.png" );
    gplotter.set_ranges( 0, 0, 0.002, 0.002 );
    gplotter.plot_png( "plasma2d_zoom.png" );

    ParticleDiagPlotter pplotter( &geom, &pdb, AXIS_X, 0.01, PARTICLE_DIAG_PLOT_HISTO2D, DIAG_R, DIAG_RP );
    pplotter.plot_png( "plasma2d_emit.png" );

    FieldDiagPlotter fplotter( &geom );
    fplotter.set_scharge( &scharge );
    fplotter.set_epot( &epot );
    fplotter.set_coordinates( 100, Vec3D(0.006,0,0), Vec3D(0.006,0.007,0) );
    field_diag_type_e diag[2] = {FIELDD_DIAG_EPOT, FIELDD_DIAG_SCHARGE};
    field_loc_type_e loc[2] = {FIELDD_LOC_Y, FIELDD_LOC_NONE};
    fplotter.set_diagnostic( diag, loc );
    fplotter.plot_png( "plasma2d_field.png" );

}


int main( int argc, char **argv )
{
    try {
	//verbose_output = 1;
	test( &argc, &argv );
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

