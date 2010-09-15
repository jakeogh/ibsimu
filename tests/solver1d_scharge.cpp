/*! \file solver1d_scharge.cpp 
 *  \brief Test solver with a 1d problem with space charge.
 *
 *  \test The simple 1d problems with constant space charge are easily
 *  solved analytically. The Poisson equation 
 *  \f[ \nabla^2 \phi = -\frac{\rho}{\epsilon} \f]
 *  can be integrated twice to get
 *  \f[ \phi = -\frac{\rho}{2 \epsilon} x^2 + Ax + B, \f]
 *  where A and B are integration constants set by boundary conditions.
 *  If we have set space charge \f$ \rho = 1\cdot10^{-4}\mathrm{~C/m}^3 \f$
 *  and Dirichlet boundaries \f$ \phi=0 \f$ at x=0 and x=0.1 m, the 
 *  constants become
 *  \f[ A=\frac{\rho}{2 \epsilon}\cdot 0.1\mathrm{~m}=564.705\mathrm{~V/m} \mathrm{~and~} B=0 \f]
 *
 */


#include <fstream>
#include <iomanip>
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"


using namespace std;


const double eps0 = 8.85418781762e-12;
double phi( double x )
{
    return( -1.0e-4*x*x/(2.0*eps0) + 1.0e-4/(2.0*eps0)*0.1*x + 0 );
}


void test( void )
{
    bool err = false;

    Geometry g( MODE_1D, Int3D(11,1,1), Vec3D(0,0,0), 0.01 );
    g.set_boundary( 1, Bound(BOUND_DIRICHLET, 0.0) );
    g.set_boundary( 2, Bound(BOUND_DIRICHLET, 0.0) );
    g.build_mesh();

    EpotProblem p;
    p.construct( g );

    ScalarField epot( g );
    ScalarField scharge( g );
    for( int a = 0; a < g.size(0); a++ )
	scharge(a) = 1.0e-4;

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    ofstream ostr( "solver1d_scharge.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "potential (V)" << "\n";
    for( int a = 0; a < g.size(0); a++ ) {
	if( fabs( epot(a) - phi(a*g.h()) ) > 1.0e-4  )
	    err = true;
	ostr << setw(14) << a*g.h() << " " 
	     << setw(14) << epot(a) << "\n";
    }
    ostr << "\n\n";
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



