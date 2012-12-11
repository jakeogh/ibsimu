/*! \file solver1d_edge.cpp 
 *  \brief Test solver with a 1d problem with smooth edges.
 *
 *  \test Test solver with a 1d problem with smooth edges.
 */

/*
 *  The simple 1d problems with constant space charge are easily
 *  solved analytically. The Poisson equation 
 *  \f[ \nabla^2 \phi = -\frac{\rho}{\epsilon} \f]
 *  can be integrated twice to get
 *  \f[ \phi = -\frac{\rho}{2 \epsilon} x^2 + Ax + B, \f]
 *  where A and B are integration constants set by boundary conditions.
 *  If we have zero space charge and Dirichlet boundaries \f$ \phi=0 \f$ at 
 *  \f$ x=x_0=0 \f$ and \f$ \phi=10\mathrm{~V} \f$ at \f$ x=x_1 \f$, 
 *  the constants become
 *  \f[ A=10\mathrm{~V} / (x_1 - x_0) \mathrm{~and~} B=0 \f]
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


double x1 = 0.013;
double x2 = 0.087;


bool s1( double x, double y, double z )
{
    return( x < x1 );
}


bool s2( double x, double y, double z )
{
    return( x > x2 );
}


double phi( double x )
{
    if( x < x1 )
	return( 0.0 );
    else if( x > x2 )
	return( 10.0 );
    return( 10.0*(x-x1)/(x2-x1) );
}


void test( int argc, char **argv )
{
    Geometry geom( MODE_1D, Int3D(11,1,1), Vec3D(0,0,0), 0.01 );
    Solid *solid1 = new FuncSolid( s1 );
    geom.set_solid( 7, solid1 );
    Solid *solid2 = new FuncSolid( s2 );
    geom.set_solid( 8, solid2 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 10.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, 10.0) );
    geom.build_mesh();

    EpotBiCGSTABSolver solver( geom );
    MeshScalarField epot( geom );
    MeshScalarField scharge( geom );

    solver.solve( epot, scharge );
    
    bool err = false;
    ofstream ostr( "solver1d_edge.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "potential (V)" << " "
	 << setw(14) << "theory (V)" << "\n";
    for( uint32_t a = 0; a < geom.size(0); a++ ) {
	double x = a*geom.h();
	if( x > x1 && x < x2 && fabs( epot(a) - phi(x) ) > 0.05  )
	    err = true;
	ostr << setw(14) << x << " " 
	     << setw(14) << epot(a) << " "
	     << setw(14) << phi(x) << "\n";
    }
    ostr << "\n\n";
    ostr.close();

    if( err ) {
	std::cout << "Error: solved potential differs from theory\n";
	exit( 1 );
    }

}


