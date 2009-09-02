/*! \file solver2d_cylindrical_beam.cpp 
 *  \brief Test solver with a 2d problem with cylindrical symmetry and beam.
 *
 *  \test The test consists of a cylindrical beam tube with a
 *  cylindrical beam inside. The test is a 2d mesh version of the
 *  problem described in solvercyl_beam.cpp.
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


bool solid1( double x, double y, double z )
{
    return( x*x+y*y >= 0.05*0.05 );
}


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
    Geometry g( MODE_2D, Int3D(121,121,1), Vec3D(0,0,0), 0.0005 );
    Solid *s1 = new FuncSolid( solid1 );
    g.set_solid( 7, s1 );
    g.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 2, Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 4, Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    g.build_mesh();

    EpotProblem p;
    p.construct( g );
    //g.debug_print();
    //p.debug_print();

    ScalarField epot( g );
    ScalarField scharge( g );
    for( int a = 0; a < g.size(0); a++ ) {
	for( int b = 0; b < g.size(1); b++ ) {
	    double x = a*g.h();
	    double y = b*g.h();
	    if( sqrt(x*x + y*y) < 0.01 )
		scharge(a,b) = 1.0e-4;
	}
    }

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    bool err = false;
    ofstream ostr( "solver2d_cylindrical_beam.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "y (m)" << " " 
	 << setw(14) << "r (m)" << " " 
	 << setw(14) << "potential (V)" << " "
	 << setw(14) << "theory (V)" << "\n";
    for( int a = 0; a < g.size(0); a++ ) {
	for( int b = 0; b < g.size(1); b++ ) {
	    double x = a*g.h();
	    double y = b*g.h();
	    double r = sqrt(x*x + y*y);
	    // Allow 10 % relative error and 10 V absolute error
	    if( r < 0.05 && fabs( (epot(a,b) - phi(r))/phi(r) ) > 0.10 &&
		fabs( (epot(a,b) - phi(r)) ) > 10.0 )
		err = true;
	    ostr << setw(14) << x << " " 
		 << setw(14) << y << " " 
		 << setw(14) << r << " " 
		 << setw(14) << epot(a,b) << " "
		 << setw(14) << phi(r) << "\n";
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

