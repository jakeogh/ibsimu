/*! \file solver1d_neumann.cpp 
 *  \brief Test solver with a 1d problem with Neumann and Dirichlet boundaries.
 *
 *  \test Test solver with a 1d problem with Neumann and Dirichlet boundaries.
 */

/*
 *  The simple 1d problems with constant space charge are easily
 *  solved analytically. The Poisson equation 
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
#include "epot_bicgstabsolver.hpp"
#include "geometry.hpp"
#include "meshscalarfield.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "ibsimu.hpp"
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


void test( int argc, char **argv )
{
    bool err1 = false;
    bool err2 = false;

    Geometry geom( MODE_1D, Int3D(11,1,1), Vec3D(0,0,0), 0.01 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,   10.0) );
    geom.build_mesh();

    MeshScalarField epot( geom );
    MeshScalarField scharge( geom );

    EpotBiCGSTABSolver solver( geom );
    solver.solve( epot, scharge );

    ofstream ostr( "solver1d_neumann.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "potential (V)" << "\n";
    for( uint32_t a = 0; a < geom.size(0); a++ ) {
	if( fabs( epot(a) - phi1(a*geom.h()) ) > 1.0e-4  )
	    err1 = true;
	ostr << setw(14) << a*geom.h() << " " 
	     << setw(14) << epot(a) << "\n";
    }
    ostr << "\n\n";

    geom.set_boundary( 1, Bound(BOUND_NEUMANN,   10.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,  0.0) );
    geom.build_mesh();
    solver.solve( epot, scharge );

    for( uint32_t a = 0; a < geom.size(0); a++ ) {
	if( fabs( epot(a) - phi2(a*geom.h()) ) > 1.0e-4  )
	    err2 = true;
	ostr << setw(14) << a*geom.h() << " " 
	     << setw(14) << epot(a) << "\n";
    }
    ostr.close();

    if( err1 ) {
	std::cout << "Error: solved potential differs from theory in case 1\n";
	exit( 1 );
    } else if( err2 ) {
	std::cout << "Error: solved potential differs from theory in case 2\n";
	exit( 1 );
    }
}



