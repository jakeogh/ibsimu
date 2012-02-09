/*! \file beam2d.cpp 
 *  \brief Test with a beam in 2d system. Space charge on axis is checked.
 *
 *  \test Test with a beam in 2d system. Space charge on axis is checked.
 */


#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "epot_gssolver.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_field.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "particledatabase.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "gtkplotter.hpp"


using namespace std;


void test( int argc, char **argv )
{    
    // 10x10 cm geometry with 0.25 cm mesh
    Geometry geom( MODE_2D, Int3D(41,41,1), Vec3D(0,0,0), 0.0025 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,  0.0) );
    geom.build_mesh();

    EpotGSSolver solver( geom );
    EpotField epot( geom );
    MeshScalarField scharge( geom );
    solver.solve( epot, scharge );

    EpotEfield efield( epot );
    MeshVectorField bfield;

    ParticleDataBase2D pdb;
    pdb.set_thread_count( 1 );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.add_2d_beam_with_energy( 100, 50.0, 1.0, 1.0, 
				 3.0e3, 0.0, 0.0, 
				 0.0, 0.0, 
				 0.0, 0.009 );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );
    solver.solve( epot, scharge );

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
    for( uint32_t j = 0; j < geom.size(1); j++ ) {

	Vec3D x( 0.05, geom.h()*j, 0.0  );
	ostr2 << setw(14) << x[1] << " "
	      << setw(14) << epot( x ) << " "
	      << setw(14) << scharge( x ) << "\n";

	for( uint32_t i = 0; i < geom.size(0); i++ ) {

	    Vec3D x( geom.h()*i, geom.h()*j, 0.0  );
	    // Check space charge density on axis
	    if( j == 0 && fabs(scharge(x)-6.57148899428e-5) > 1e-8 )
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
    ostr2.close();

    if( false ) {
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.set_scharge( &scharge );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
    }

    if( err ) {
	std::cout << "Error: calculated space charge differs from theory\n";
	exit( 1 );
    }
}
