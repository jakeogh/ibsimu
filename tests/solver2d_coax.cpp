/*! \file solver2d_coax.cpp 
 *  \brief Test solver with a 2d problem made of two concentric cylinders.
 *
 *  \test Test solver with a 2d problem made of two concentric cylinders.
 */

/*
 *  Analytically this problem is solved by using cylindrical coordinates.
 *  The Poisson equation 
 *  \f[ \nabla^2 \phi = -\frac{\rho}{\epsilon} \f]
 *  becomes
 *  \f[ \frac{\partial^2 \phi}{\partial r^2} + 
    \frac{1}{r} \frac{\partial \phi}{\partial r} +
    \frac{1}{r^2} \frac{\partial^2 \phi}{\partial \theta^2} +
    \frac{\partial^2 \phi}{\partial z^2} = -\frac{\rho}{\epsilon}. \f]
 *  The third and fourth terms are zero in this symmetric case and in case
 *  of no space charge, the differential equation becomes
 *  \f[ \frac{\partial^2 \phi}{\partial r^2} + 
    \frac{1}{r} \frac{\partial \phi}{\partial r} = 0. \f]
 *  By setting \f$ \frac{\partial \phi}{\partial r} = y \f$ and integrating
 *  the resulting separable first order DE we get  
 *  \f[ \frac{\partial \phi}{\partial r} = A \frac{1}{r}, \f]
 *  which can again be integrated to get the final solution
 *  \f[ \phi = A \ln r + B, \f]
 *  where A and B are integration constants set by boundary conditions.
 *  If we have Dirichlet boundaries \f$ \phi=0\mathrm{~V} \f$ at 
 *  \f$ x=0.02\mathrm{~m} \f$ and \f$ \phi=10\mathrm{~V} \f$ at 
 *  \f$ x=0.07\mathrm{~m} \f$, the constants become
 *  \f[ A=7.98235600148\mathrm{~V} \mathrm{~and~} B=31.2271603153 \f]
 *
 */


#include <fstream>
#include <iomanip>
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "gtkplotter.hpp"
#include "error.hpp"
#include "ibsimu.hpp"


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

    EpotProblem p;
    p.construct( geom );
    //geom.debug_print();
    //p.debug_print();

    ScalarField epot( geom );
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    bool err = false;
    ofstream ostr( "solver2d_coax.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "y (m)" << " " 
	 << setw(14) << "r (m)" << " " 
	 << setw(14) << "potential (V)" << " "
	 << setw(14) << "theory (V)" << "\n";
    for( int a = 0; a < geom.size(0); a++ ) {
	for( int b = 0; b < geom.size(1); b++ ) {
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

    /*
    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.new_geometry_plot_window();
    plotter.run();
    */

    if( err ) {
	std::cout << "Error: solved potential differs from theory\n";
	exit( 1 );
    }
}


