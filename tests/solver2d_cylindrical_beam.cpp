/*! \file solver2d_cylindrical_beam.cpp 
 *  \brief Test solver with a 2d problem with cylindrical symmetry and beam.
 *
 *  \test Test solver with a 2d problem with cylindrical symmetry and beam.
 */


#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "geomplotter.hpp"
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"


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


void test( int argc, char **argv )
{
    Geometry geom( MODE_2D, Int3D(121,121,1), Vec3D(0,0,0), 0.0005 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    geom.build_mesh();

    EpotProblem p;
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );
    for( int a = 0; a < geom.size(0); a++ ) {
	for( int b = 0; b < geom.size(1); b++ ) {
	    double x = a*geom.h();
	    double y = b*geom.h();
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
    for( int a = 0; a < geom.size(0); a++ ) {
	for( int b = 0; b < geom.size(1); b++ ) {
	    double x = a*geom.h();
	    double y = b*geom.h();
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


