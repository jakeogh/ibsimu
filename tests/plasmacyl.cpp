/*! \file plasmacyl.cpp
 *  \brief Test with a plasma in cylindrical electrode configuration.
 *
 *  \test Test with a plasma in cylindrical electrode configuration.
 */


#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "particledatabase.hpp"
#include "geometry.hpp"
#include "convergence.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "vectorfield.hpp"
#include "meshvectorfield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "particlediagplotter.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"


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


bool init_plasma( double x, double y, double z )
{
    return( x < 0.00055 );
}


void test( int argc, char **argv )
{
    // 12x7 mm geometry with 0.05 mm mesh size
    Geometry geom( MODE_CYL, Int3D(241,141,1), Vec3D(0,0,0), 0.00005 );
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
    p.set_initial_plasma( 5.0, init_plasma );
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );
    ScalarField scharge_ave( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );

    MeshVectorField bfield;
    EpotEfield efield( geom, epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_MIRROR,      FIELD_EXTRAPOLATE,
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBaseCyl pdb;
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );

    Convergence conv;
    conv.add_epot( epot, 1, 1, 1.0e-6 );
    conv.add_scharge( scharge_ave, 1, 1, 1.0e-6 );
    conv.add_tdiag( pdb, AXIS_X, 11.9e-3, 1, 1, 1.0e-6 );

    for( size_t i = 0; i < 2; i++ ) {

	if( i == 1 ) {
	    double rhoe = pdb.get_rhosum();
	    p.set_pexp_plasma( -rhoe, 5.0, 5.0 );
	    p.construct( geom );
	}

	p.solve( epot, scharge_ave );

	pdb.clear();
	pdb.add_2d_beam_with_energy( 5000, 600.0, 1.0, 1.0, 
				     5.0, 0.0, 0.5, 
				     0.0, 0.0, 
				     0.0, 0.0015 );
	pdb.iterate_trajectories( scharge, efield, bfield, geom );

	if( i == 0 ) {
            scharge_ave = scharge;
        } else {
            double coef = 0.3;
            scharge *= coef;
            scharge_ave += scharge;
            scharge_ave *= (1.0/(1.0+coef));
        }

	conv.evaluate_iteration();

	if( false ) {
	    ScalarField tdens( geom );
	    pdb.build_trajectory_density_field( tdens );
	    GTKPlotter plotter( &argc, &argv );
	    plotter.set_geometry( &geom );
	    plotter.set_epot( &epot );
	    plotter.set_trajdens( &tdens );
	    plotter.set_scharge( &scharge );
	    plotter.set_particledatabase( &pdb );
	    plotter.new_geometry_plot_window();
	    plotter.run();
	}
    }

    ofstream ofconv( "plasmacyl_conv.dat" );
    conv.print_history( ofconv );
    ofconv.close();

    ScalarField tdens( geom );
    pdb.build_trajectory_density_field( tdens );

    GeomPlotter gplotter( &geom );
    gplotter.set_size( 1024, 768 );
    gplotter.set_epot( &epot );
    std::vector<double> eqlines;
    eqlines.push_back( -8.0 );
    eqlines.push_back( -4.0 );
    eqlines.push_back( 0.01 );
    eqlines.push_back( +4.0 );
    gplotter.set_eqlines_manual( eqlines );
    gplotter.set_particle_database( &pdb );
    gplotter.set_particle_div( 0 );
    gplotter.set_trajdens( &tdens );
    gplotter.set_fieldgraph_plot( FIELD_TRAJDENS );
    gplotter.plot_png( "plasmacyl.png" );

    ParticleDiagPlotter pplotter1( &geom, &pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_R, DIAG_RP );
    pplotter1.plot_png( "plasmacyl1_r_rp.png" );

    ParticleDiagPlotter pplotter2( &geom, &pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_R, DIAG_AP );
    pplotter2.plot_png( "plasmacyl1_r_ap.png" );

    ParticleDiagPlotter pplotter3( &geom, &pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_RP, DIAG_AP );
    pplotter3.plot_png( "plasmacyl1_rp_ap.png" );



    ParticleDiagPlotter pplotter1b( &geom, &pdb, AXIS_X, 0.012-1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_R, DIAG_RP );
    pplotter1b.plot_png( "plasmacyl2_r_rp.png" );

    ParticleDiagPlotter pplotter2b( &geom, &pdb, AXIS_X, 0.012-1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_R, DIAG_AP );
    pplotter2b.plot_png( "plasmacyl2_r_ap.png" );

    ParticleDiagPlotter pplotter3b( &geom, &pdb, AXIS_X, 0.012-1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_RP, DIAG_AP );
    pplotter3b.plot_png( "plasmacyl2_rp_ap.png" );
}

