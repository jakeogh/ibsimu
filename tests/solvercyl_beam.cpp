/*! \file solvercyl_beam.cpp 
 *  \brief Test solver with a cylindrical problem with localized space charge.
 *
 *  \test Test solver with a cylindrical problem with localized space charge.
 */


#include <fstream>
#include <iomanip>
#include "epot_bicgstabsolver.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "constants.hpp"


using namespace std;


double r1 = 0.01;
double r2 = 1000.0*5.0e-5;
double rho = 1.0e-4;
double I_over_v = M_PI*r1*r1*rho;
double coef = I_over_v/(2.0*M_PI*EPSILON0);



double phi( double r )
{
    if( r < 0.01 )
	return( -coef*( r*r/(2.0*r1*r1) + log(r1/r2) - 0.5 ) );
    else if( r < 0.05 )
	return( -coef*log(r/r2) );
    else
	return( 0.0 );
}


void test( int argc, char **argv )
{
    bool err = false;

    Geometry geom( MODE_CYL, Int3D(5,1001,1), Vec3D(0,0,0), 5.0e-5 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,   0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,   0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,   0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET, 0.0) );
    geom.build_mesh();

    EpotBiCGSTABSolver solver( geom );
    EpotField epot( geom );
    MeshScalarField scharge( geom );
    for( unsigned int a = 0; a < geom.size(0); a++ ) {
	for( unsigned int b = 0; b < geom.size(1); b++ ) {
	    //double x = a*geom.h();
	    double r = b*geom.h();
	    if( r < 0.01 )
		scharge(a,b) = 1.0e-4;
	}
    }

    solver.solve( epot, scharge );

    ofstream ostr( "solvercyl_beam.dat" );
    ostr << "# "
	 << setw(12) << "r (m)" << " " 
	 << setw(14) << "potential (V)" << " "
	 << setw(14) << "theory (V)" << "\n";
    for( unsigned int a = 0; a < geom.size(0); a++ ) {
	for( unsigned int b = 0; b < geom.size(1); b++ ) {
	    //double x = a*geom.h();
	    double r = b*geom.h();
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



