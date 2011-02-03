/*! \file solvercyl_beam.cpp 
 *  \brief Test solver with a cylindrical problem with localized space charge.
 *
 *  \test Test solver with a cylindrical problem with localized space charge.
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


void test( int argc, char **argv )
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



