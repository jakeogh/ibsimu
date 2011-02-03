/*! \file gs_solver.cpp 
 *  \brief Test Gauss-Seidel solver with a 2d problem made of two concentric cylinders.
 *
 *  \test Test Gauss-Seidel solver with a 2d problem made of two concentric cylinders.
 */

#include <fstream>
#include <iomanip>
#include "gs_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
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


void check_epot( Geometry &geom, ScalarField &epot )
{
    bool err = false;
    ofstream ostr( "gs_solver.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "y (m)" << " " 
	 << setw(14) << "r (m)" << " " 
	 << setw(14) << "potential (V)" << " "
	 << setw(14) << "theory (V)" << "\n";
    double maxerr = 0.0;
    int maxerrl[2] = {0,0};
    for( int a = 0; a < geom.size(0); a++ ) {
	for( int b = 0; b < geom.size(1); b++ ) {
	    double x = a*geom.h();
	    double y = b*geom.h();
	    double r = sqrt(x*x + y*y);
	    double error = fabs( epot(a,b) - phi(r) );
	    if( r > 0.02 && r < 0.07 && error > 0.5  )
		err = true;
	    if( r > 0.02 && r < 0.07 && error > maxerr  ) {
		maxerr = error;
		maxerrl[0] = a;
		maxerrl[1] = b;
	    }
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
	std::cout << "Maximum error = " << maxerr << "\n";
	std::cout << " at (" 
		  << maxerrl[0] << ", " 
		  << maxerrl[1] << ")\n";
	exit( 1 );
    }
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
    p.set_neumann_order( 1 );
    p.enable_smooth_solids( false );
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );
    GSSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );
    check_epot( geom, epot );
}

