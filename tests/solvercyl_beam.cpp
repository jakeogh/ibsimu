/*! \file solvercyl_beam.cpp 
 *  \brief Test solver with a cylindrical problem with localized space charge.
 *
 *  \test The test consists of a cylindrical beam tube with radius 
 *  \f$ r_1 = 0.05\mathrm{~cm} \f$. There is a beam with radius \f$ r_0 = 
    0.01\mathrm{~cm} \f$ in the tube with constant space charge 
 *  \f$ \rho = 1\cdot10^{-4}\mathrm{~C/m}^3 \f$ inside. The problem can
 *  be solved analytically by integrating the Poisson equation in the 
 *  two domains and by using the boundary conditions. The solution is
 *  \f{eqnarray*} \phi &=& -\frac{\rho}{4 \epsilon_0}r^2 + C, r <= r_0 \\
                  \phi &=& A \log r + B, r > r_0, \f}
 *  where the constants are 
 *  \f[ A = -\frac{\rho r_0^2}{2 \epsilon_0} \mathrm{~and~} 
        B = -A \log r_1 \mathrm{~and~} C = A \log \frac{r_0}{r_1} + 
        \frac{\rho r_0^2}{4 \epsilon_0} \f]
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
double phi( double r )
{
    if( r < 0.01 )
	return( -2823522.6669*r*r + 1191.20915203 );
    else if( r < 0.05 )
	return( -564.704533379*log(r) - 1691.70359567 );
    else
	return( 0.0 );
}


void test( void )
{
    bool err = false;

    Geometry g( MODE_CYL, Int3D(5,1001,1), Vec3D(0,0,0), 5.0e-5 );
    g.set_boundary( 1, Bound(BOUND_NEUMANN,   0.0) );
    g.set_boundary( 2, Bound(BOUND_NEUMANN,   0.0) );
    g.set_boundary( 3, Bound(BOUND_NEUMANN,   0.0) );
    g.set_boundary( 4, Bound(BOUND_DIRICHLET, 0.0) );
    g.build_mesh();

    EpotProblem p;
    p.construct( g );

    ScalarField epot( g );
    ScalarField scharge( g );
    for( int a = 0; a < g.size(0); a++ ) {
	for( int b = 0; b < g.size(1); b++ ) {
	    //double x = a*g.h();
	    double r = b*g.h();
	    if( r < 0.01 )
		scharge(a,b) = 1.0e-4;
	}
    }

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    ofstream ostr( "solvercyl_beam.dat" );
    ostr << "# "
	 << setw(12) << "r (m)" << " " 
	 << setw(14) << "potential (V)" << " "
	 << setw(14) << "theory (V)" << "\n";
    for( int a = 0; a < g.size(0); a++ ) {
	for( int b = 0; b < g.size(1); b++ ) {
	    //double x = a*g.h();
	    double r = b*g.h();
	    if( r < 0.05 && fabs( (epot(a,b) - phi(r))/phi(r) ) > 0.05 &&
		fabs( (epot(a,b) - phi(r)) ) > 0.1 )
		err = true;
	    ostr << setw(14) << r << " " 
		 << setw(14) << epot(a,b) << " "
		 << setw(14) << phi(r) << "\n";
	}
	ostr << "\n\n";
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

