/*! \file solver1d_hard_edge.cpp 
 *  \test  Test solver with a 1d problem with hard edges.
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


bool s1( double x, double y, double z )
{
    return( x > 0.095 );
}


bool s2( double x, double y, double z )
{
    return( x > 0.085 );
}


double phi1( double x )
{
    return( 10.0/0.1*x );
}


double phi2( double x )
{
    if( x < 0.09 )
	return( 10.0/0.09*x );
    else
	return( 10.0 );
}


void test( void )
{
    bool err = false;

    Geometry g( MODE_1D, Int3D(11,1,1), Vec3D(0,0,0), 0.01 );
    Solid *solid1 = new FuncSolid( s1 );
    g.set_solid( 7, solid1 );
    g.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 2, Bound(BOUND_DIRICHLET, 10.0) );
    g.set_boundary( 7, Bound(BOUND_DIRICHLET, 10.0) );
    g.build_mesh();

    EpotProblem p;
    p.enable_smooth_solids( false );
    p.construct( g );

    ScalarField epot( g );
    ScalarField scharge( g );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );
    
    ofstream ostr( "solver1d_hard_edge1.dat" );
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

    Solid *solid2 = new FuncSolid( s2 );
    g.set_solid( 7, solid2 );
    g.build_mesh();
    p.construct( g );
    p.solve( epot, scharge );

    ostr.open( "solver1d_hard_edge2.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "potential (V)" << "\n";
    for( int a = 0; a < g.size(0); a++ ) {
	if( fabs( epot(a) - phi2(a*g.h()) ) > 1.0e-4  )
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

