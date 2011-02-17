/*! \file solver3d_sphere.cpp 
 *  \brief Test solver with a 3d problem made of two concentric spheres.
 *
 *  \test Test solver with a 3d problem made of two concentric spheres.
 */


#include <fstream>
#include <iomanip>
#include "epot_gssolver.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"

#include "gtkplotter.hpp"

using namespace std;


double phi_a = 0.0;
double phi_b = 10.0;
double r_a = 0.02;
double r_b = 0.07;
double A = (phi_a - phi_b)/(1.0/r_a - 1.0/r_b);
double B = phi_a + A/r_a;


bool solid1( double x, double y, double z )
{
    return( x*x+y*y+z*z <= 0.02*0.02 );
}


bool solid2( double x, double y, double z )
{
    return( x*x+y*y+z*z >= 0.07*0.07 );
}


double phi( double r )
{
    if( r <= 0.02 )
	return( phi_a );
    else if( r >= 0.07 )
	return( phi_b );
    return( A/r - B );
}

/* Optimize w at h=0.001
 *
 * w = 1.00, iter = 7380
 *
 * w = 1.85, iter = 755
 * w = 1.90, iter = 484
 * w = 1.91, iter = 424
 * w = 1.92, iter = 359
 * w = 1.93, iter = 279
 * w = 1.94, iter = 315
 * w = 1.95, iter = 375
 * w = 1.96, iter = 465
 * w = 1.99, iter = 1794
 *
 */

/* w = 1.93, -O2 compiler
 *
 * h = 0.002,   iter = 260, time = 0.26 s
 * h = 0.001,   iter = 279, time = 2.01 s
 * h = 0.0005,  iter = 1272, time = 68.58 s
 * h = 0.00025, iter = 4542, time = 1993.02 s
 *
 *
 *
 */

void test( int argc, char **argv )
{
    double h = 0.002;
    int32_t size = (int32_t)ceil(0.08/h) + 1;
    Geometry g( MODE_3D, Int3D(size,size,size), Vec3D(0,0,0), h );
    Solid *s1 = new FuncSolid( solid1 );
    g.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    g.set_solid( 8, s2 );
    g.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 2, Bound(BOUND_DIRICHLET, 10.0) );
    g.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 4, Bound(BOUND_DIRICHLET, 10.0) );
    g.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 8, Bound(BOUND_DIRICHLET, 10.0) );
    g.build_mesh();

    EpotGSSolver solver( g );
    MeshScalarField epot( g );
    MeshScalarField scharge( g );

    solver.set_imax( 10000 );
    solver.set_w( 1.93 );
    solver.solve( epot, scharge );

    bool err = false;
    double maxerr = 0.0;
    int maxerrl[3] = {0,0,0};
    ofstream ostr( "solver3d_sphere.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "y (m)" << " " 
	 << setw(14) << "z (m)" << " " 
	 << setw(14) << "r (m)" << " " 
	 << setw(14) << "potential (V)" << " "
	 << setw(14) << "theory (V)" << "\n";
    for( uint32_t a = 0; a < g.size(0); a++ ) {
	for( uint32_t b = 0; b < g.size(1); b++ ) {
	    for( uint32_t c = 0; c < g.size(2); c++ ) {
		double x = a*g.h();
		double y = b*g.h();
		double z = c*g.h();
		double r = sqrt(x*x + y*y + z*z);
		double e = epot(a,b,c);
		double error = fabs( epot(a,b,c) - phi(r) );
		if( r < 0.02 )
		    e = 0.0;
		else if( r > 0.07 )
		    e = 10.0;
		else if( r > 0.025 && error > 0.32 )
		    err = true;
		if( r > 0.02 && r < 0.07 && error > maxerr ) {
		    maxerr = error;
		    maxerrl[0] = a;
		    maxerrl[1] = b;
		    maxerrl[2] = c;
		}
		ostr << setw(14) << x << " " 
		     << setw(14) << y << " " 
		     << setw(14) << z << " " 
		     << setw(14) << r << " " 
		     << setw(14) << e << " "
		     << setw(14) << phi(r) << "\n";
	    }
	}
    }

    ostr.close();

    /*
    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &g );
    plotter.set_scharge( &scharge );
    plotter.set_epot( &epot );
    plotter.new_geometry_plot_window();
    plotter.run();
    */

    if( err ) {
	/*
	  Constructing geometry
	  origo =            0            0            0
	  size  =           41           41           41
	  max   =         0.08         0.08         0.08
	  h     = 0.002
	  Constructing linear electric potential problem
	  dof = 25098
	  Solving problem
	  Using UMFPACK solver
	  Done
	  time used = 3.88 s (3.88477 s realtime)
	  Maximum error = 0.305325
	  at (0, 0, 13)
	*/
	std::cout << "Error: solved potential differs from theory\n";
	std::cout << "Maximum error = " << maxerr << "\n";
	std::cout << " at (" 
		  << maxerrl[0] << ", " 
		  << maxerrl[1] << ", " 
		  << maxerrl[2] << ")\n";
	exit( 1 );
    }
}


