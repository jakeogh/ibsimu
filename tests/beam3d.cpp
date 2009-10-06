/*! \file beam3d.cpp 
 *  \test Test with a beam in 3d system.
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


void test( int *argc, char ***argv )
{
    // 10x20x40 cm geometry with 0.25 cm mesh
    Geometry geom( MODE_3D, Int3D(41,41,41), Vec3D(0,0,0), 0.0025 );
    //Geometry geom( MODE_3D, Int3D(51,101,101), Vec3D(0,-0.1,-0.1), 0.002 );

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
    VectorField bfield;

    ParticleDataBase3D pdb;
    pdb.set_thread_count( 4 );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );
    for( size_t a = 0; a < 1; a++ ) {

	p.solve( epot, scharge );

	pdb.clear();
	pdb.add_cylindrical_beam_with_energy( 1000, 50.0, 1.0, 1.0, 
					      3.0e3, 0.0, 0.5,
					      Vec3D(0,0,0), // center
					      Vec3D(0,1,0), // dir1
					      Vec3D(0,0,1), // dir2
					      0.005 );      // radius
	pdb.iterate_trajectories( scharge, efield, bfield, geom );
    }

    /*
    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_scharge( &scharge );
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
    gplotter.set_view( VIEW_XY, 0 );
    gplotter.plot_png( "beam3d_xy.png" );
    gplotter.set_view( VIEW_YZ, 0 );
    gplotter.plot_png( "beam3d_yz.png" );
}


int main( int argc, char **argv )
{
    try {
	verbose_output = 0;
    	test( &argc, &argv );
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

