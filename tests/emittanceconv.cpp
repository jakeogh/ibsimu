/*! \file emittanceconv.cpp 
 *  \test Test emittance (r,r') to (x,x') converter
 */


#include <fstream>
#include <iomanip>
#include <iostream>
#include "geometry.hpp"
#include "epot_problem.hpp"
#include "epot_efield.hpp"
#include "bicgstab_solver.hpp"
#include "particledatabase.hpp"
#include "trajectorydiagnostics.hpp"
#include "verbose.hpp"
#include "error.hpp"

#include "frame.hpp"
#include "colormap.hpp"
#include "palette.hpp"


using namespace std;


void test( void )
{
    Geometry geom( MODE_CYL, Int3D(11,11,1), Vec3D(0,0,0), 1e-3 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 3, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET, 0.0) );
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

    // Make a circular beam of Neon (m=20) 6+ with J = 45 A/m2 and radius 5 mm.
    // Starting energy is 3 keV and transverse temperature 1.0 eV.
    // The emittance should is known analytically and is 0.5*sqrt(kT/m)*r/v_z = 3.22748612e-5 pi*m*rad 
    ParticleDataBaseCyl pdb;
    pdb.add_2d_beam_with_energy( 20000, 45.0, 6.0, 20.0, 
				 3.0e3, 0.0, 1.0, 
				 0.0, 0.0, 
				 0.0, 0.005 );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    // Get particle diagnostics from beam
    TrajectoryDiagnosticData tdata;
    vector<trajectory_diagnostic_e> diagnostics;
    diagnostics.push_back( DIAG_R );
    diagnostics.push_back( DIAG_RP );
    diagnostics.push_back( DIAG_AP );
    diagnostics.push_back( DIAG_CURR );
    pdb.trajectories_at_plane( tdata, AXIS_X, 1e-6, diagnostics );

    // Make emittance conversion
    EmittanceConv emit( 110, 110, tdata(0).data(), tdata(1).data(), tdata(2).data(), tdata(3).data() );
    std::cout << "xave    = " << emit.xave() << "\n";
    std::cout << "xpave   = " << emit.xpave() << "\n";
    std::cout << "alpha   = " << emit.alpha() << "\n";
    std::cout << "beta    = " << emit.beta() << "\n";
    std::cout << "gamma   = " << emit.gamma() << "\n";
    std::cout << "epsilon = " << emit.epsilon() << "\n";

    // Get emittance histogram and plot it
    const Histogram2D &histo = emit.histogram();

    // Make plot surface and frame
    cairo_surface_t *surface;
    cairo_t *cairo;
    Frame frame;

    // Colormap
    double range[4];
    histo.get_range( range );
    Colormap colormap( range, histo.n(), histo.m(), histo.get_data() );
    Palette palette;
    palette.clear();
    palette.push_back( Color(1,1,1), 0 );
    palette.push_back( Color(1,1,0), 1 );
    palette.push_back( Color(1,0,0), 2 );
    palette.push_back( Color(0,0,0), 3 );
    palette.norm();
    colormap.set_palette( palette );
    colormap.set_interpolation( INTERPOLATION_CLOSEST );
    frame.add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, &colormap );
    
    // Plot frame and colormap
    surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, 640, 480 );
    cairo = cairo_create( surface );
    frame.draw( cairo );
    cairo_surface_write_to_png( surface, "emittanceconv.png" );
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );

    // Output emittance data
    ofstream ostr;
    ostr.open( "emittanceconv.dat" );
    for( int j = histo.m()-1; j >= 0; j-- ) {
	for( int i = 0; i < (int)histo.n(); i++ )
	    ostr << setw(12) << histo(i,j) << " ";
	ostr << "\n";
    }
    ostr.close();
}


int main( void )
{
    try {
	verbose_output = 1;
	test();
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

