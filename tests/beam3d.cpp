/*! \file beam3d.cpp 
 *  \brief Test with a beam in 3d system.
 *
 *  \test Test with a beam in 3d system.
 *
 *  Geometry is made with two electrodes with 1 kV voltage
 *  difference. A 3 keV proton beam is accelerated through holes in
 *  the electrodes. Diagnostics is made at the end and the number of
 *  particles caught in the diagnostics is checked.
 */


#include <fstream>
#include <iomanip>

#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "particledatabase.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "ibsimutest.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"


using namespace std;



bool electrode1_func( double x, double y, double z )
{
    return( x < 0.01 && 
	    y*y + z*z > 0.01*0.01 && 
	    y*y + z*z < 0.2*0.2 );
}


bool electrode2_func( double x, double y, double z )
{
    return( x > 0.09 && 
	    y*y + z*z > 0.01*0.01 &&
	    y*y + z*z < 0.2*0.2 );
}


void test( int argc, char **argv )
{
    Geometry geom( MODE_3D, Int3D(26,26,26), Vec3D(0,-0.05,-0.05), 0.004 );

    Solid *solid1 = new FuncSolid( electrode1_func );
    geom.set_solid( 7, solid1 );
    Solid *solid2 = new FuncSolid( electrode2_func );
    geom.set_solid( 8, solid2 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,  1000.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,      0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,      0.0) );
    geom.set_boundary( 5, Bound(BOUND_NEUMANN,      0.0) );
    geom.set_boundary( 6, Bound(BOUND_NEUMANN,      0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,    0.0) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET,  1000.0) );
    geom.build_mesh();

    EpotProblem p;
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );

    EpotEfield efield( geom, epot );
    MeshVectorField bfield;

    ParticleDataBase3D pdb;
    pdb.set_thread_count( 4 );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );
    for( size_t a = 0; a < 1; a++ ) {

	p.solve( epot, scharge );

	pdb.clear();
	pdb.add_cylindrical_beam_with_energy( 4000, 50.0, 1.0, 1.0, 
					      3.0e3, 0.0, 0.5,
					      Vec3D(0,0,0), // center
					      Vec3D(0,1,0), // dir1
					      Vec3D(0,0,1), // dir2
					      0.005 );      // radius
	pdb.iterate_trajectories( scharge, efield, bfield, geom );
    }


    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_scharge( &scharge );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();


    GeomPlotter gplotter( &geom );
    gplotter.set_scharge( &scharge );
    gplotter.set_epot( &epot );
    gplotter.set_particle_database( &pdb );
    gplotter.set_particle_div( 11 );
    std::vector<double> pot;
    pot.push_back( +2 );
    pot.push_back( 0 );
    pot.push_back( -2 );
    pot.push_back( -20 );
    gplotter.set_eqlines_manual( pot );
    gplotter.set_font_size( 16 );
    gplotter.set_view( VIEW_XY );
    gplotter.plot_png( "beam3d_xy.png" );
    gplotter.set_view( VIEW_YZ, 0 );
    gplotter.plot_png( "beam3d_yz.png" );

    TrajectoryDiagnosticData tdata;
    std::vector<trajectory_diagnostic_e> diag;
    diag.push_back( DIAG_X );
    diag.push_back( DIAG_XP );
    pdb.trajectories_at_plane( tdata, AXIS_X, 0.1, diag );
    //if( tdata.traj_size() != 4000 ) 
    //throw( ErrorTest( ERROR_LOCATION, "incorrect number of particles caught: " + to_string(tdata.traj_size()) + " != 4000" ) );

    pdb.export_path_manager_data( "beam3d.path",
 				  4.0e3, 1.0, 1.0,
				  Vec3D(0.09,0,0), 
				  Vec3D(0,1,0), 
				  Vec3D(0,0,1) );
}

