/*! \file solver3d_sphere.cpp 
 *  \brief Test solver with a 3d problem made of two concentric spheres.
 *
 *  \test Analytically this problem is solved by using spherical coordinates.
 *  The Poisson equation 
 *  \f[ \nabla^2 \phi = -\frac{\rho}{\epsilon} \f]
 *  becomes
 *  \f[ \frac{1}{r^2} \frac{\partial}{\partial r}
    \left( r^2 \frac{\partial \phi}{\partial r} \right)
    = -\frac{\rho}{\epsilon}. \f]
 *  In caseof no space charge, the solution is
 *  \f[ \phi = \frac{A}{r} - B, \f]
 *  where \f$ A = \frac{\phi_a - \phi_b}{\frac{1}{r_a} - \frac{1}{r_b}} \f$ and 
 *  \f$ B = \phi_a + \frac{A}{r_a} \f$. 
 *
 */


#include <fstream>
#include <iomanip>
#include "fileplot.hpp"
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"


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


void test( void )
{
    verbose_output = 1;

    Geometry g( MODE_3D, Int3D(41,41,41), Vec3D(0,0,0), 0.002 );
    Solid *s1 = new FuncSolid( solid1 );
    g.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    g.set_solid( 8, s2 );
    g.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 2, Bound(BOUND_DIRICHLET, 10.0) );
    g.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 4, Bound(BOUND_DIRICHLET, 10.0) );
    g.set_boundary( 7, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 8, Bound(BOUND_DIRICHLET, 10.0) );
    g.build_mesh();

    EpotProblem p;
    p.construct( g );
    //g.debug_print();
    //p.debug_print();

    ScalarField epot( g );
    ScalarField scharge( g );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    bool err = false;
    ofstream ostr( "solver3d_sphere.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "y (m)" << " " 
	 << setw(14) << "z (m)" << " " 
	 << setw(14) << "r (m)" << " " 
	 << setw(14) << "potential (V)" << " "
	 << setw(14) << "theory (V)" << "\n";
    for( int a = 0; a < g.size(0); a++ ) {
	for( int b = 0; b < g.size(1); b++ ) {
	    for( int c = 0; c < g.size(2); c++ ) {
		double x = a*g.h();
		double y = b*g.h();
		double z = c*g.h();
		double r = sqrt(x*x + y*y + z*z);
		if( r > 0.02 && r < 0.07 && fabs( epot(a,b,c) - phi(r) ) > 0.03  )
		    err = true;
		ostr << setw(14) << x << " " 
		     << setw(14) << y << " " 
		     << setw(14) << z << " " 
		     << setw(14) << r << " " 
		     << setw(14) << epot(a,b,c) << " "
		     << setw(14) << phi(r) << "\n";
	    }
	}
    }

    ostr.close();

    if( err ) {
	std::cout << "Error: solved potential differs from theory\n";
	exit( 1 );
    }
}


int main( void )
{
    try {
	test();
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

