/*! \file solver1d_neumann.cpp 
 *  \brief Test solver with a 1d problem with Neumann and Dirichlet boundaries.
 *
 *  \test The simple 1d problems with constant space charge are easily
 *  solved analyticalle. The Poisson equation 
 *  \f[ \nabla^2 \phi = -\frac{\rho}{\epsilon} \f]
 *  can be integrated twice to get
 *  \f[ \phi = -\frac{\rho}{2 \epsilon} x^2 + Ax + B, \f]
 *  where A and B are integration constants set by boundary conditions.
 *  If we have zero space charge and Dirichlet boundary \f$ \phi=0 \f$ at 
 *  x=0 and Neumann boundary \f$ d\phi/dn = 10~\mathrm{V/m} \f$ at x=0.1 m, 
 *  the constants become
 *  \f[ A=-10\mathrm{~V/m} \mathrm{~and~} B=0 \f]
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
double phi1( double x )
{
    return( -10.0*x + 0.0 );
}
double phi2( double x )
{
    return( 10.0*x - 1.0 );
}


void test( void )
{
    bool err = false;

    Geometry g( MODE_1D, Int3D(11,1,1), Vec3D(0,0,0), 0.01 );
    g.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 2, Bound(BOUND_NEUMANN,   10.0) );
    g.build_mesh();

    EpotProblem p;
    p.construct( g );

    ScalarField epot( g );
    ScalarField scharge( g );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    ofstream ostr( "solver1d_neumann.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "potential (V)" << "\n";
    for( int a = 0; a < g.size(0); a++ ) {
	if( fabs( epot(a) - phi1(a*g.h()) ) > 1.0e-4  )
	    err = true;
	ostr << setw(14) << a*g.h() << " " 
	     << setw(14) << epot(a) << "\n";
    }
    ostr << "\n\n";
    ostr.close();

    g.set_boundary( 1, Bound(BOUND_NEUMANN,   10.0) );
    g.set_boundary( 2, Bound(BOUND_DIRICHLET,  0.0) );
    g.build_mesh();
    p.construct( g );
    p.solve( epot, scharge );

    for( int a = 0; a < g.size(0); a++ ) {
	if( fabs( epot(a) - phi2(a*g.h()) ) > 1.0e-4  )
	    err = true;
    }

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

