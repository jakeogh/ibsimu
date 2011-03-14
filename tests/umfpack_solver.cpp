/*! \file tests/umfpack_solver.cpp 
 *  \brief Test UMFPACK solver with a 2d problem made of two concentric cylinders.
 *
 *  \test Test UMFPACK solver with a 2d problem made of two concentric cylinders.
 */

#include <fstream>
#include <iomanip>
#include "epot_umfpacksolver.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "gtkplotter.hpp"
#include "ibsimu.hpp"
#include "ibsimutest.hpp"


using namespace std;


bool solid1( double x, double y, double z )
{
    return( x*x+y*y <= 0.02*0.02 );
}


bool solid2( double x, double y, double z )
{
    return( x*x+y*y >= 0.07*0.07 );
}


double phi( double r )
{
    return( 7.98235600148*log(r) + 31.2271603153 );
}


void test( int argc, char **argv )
{
    Geometry geom( MODE_2D, Int3D(41,41,1), Vec3D(0,0,0), 0.002 );
    //Geometry geom( MODE_2D, Int3D(6,6,1), Vec3D(0,0,0), 0.016 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 10.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET, 10.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, 10.0) );
    geom.build_mesh();

    EpotUMFPACKSolver solver( geom );
    EpotField epot( geom );
    MeshScalarField scharge( geom );
    solver.solve( epot, scharge );

    bool err = false;
    ofstream ostr( "umfpack_solver.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "y (m)" << " " 
	 << setw(14) << "r (m)" << " " 
	 << setw(14) << "potential (V)" << " "
	 << setw(14) << "theory (V)" << "\n";
    for( uint32_t a = 0; a < geom.size(0); a++ ) {
	for( uint32_t b = 0; b < geom.size(1); b++ ) {
	    double x = a*geom.h();
	    double y = b*geom.h();
	    double r = sqrt(x*x + y*y);
	    if( r > 0.02 && r < 0.07 && fabs( epot(a,b) - phi(r) ) > 0.15  )
		err = true;
	    ostr << setw(14) << x << " " 
		 << setw(14) << y << " " 
		 << setw(14) << r << " " 
		 << setw(14) << epot(a,b) << " "
		 << setw(14) << phi(r) << "\n";
	}
    }

    ostr.close();

    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.new_geometry_plot_window();
    plotter.run();

    if( err )
	throw( ErrorTest( ERROR_LOCATION, "Error: solved potential differs from theory" ) );
}

