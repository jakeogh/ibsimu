/*! \file efield.cpp 
 *  \brief Test for electric field evaluator in 1d.
 *
 *  \test Test for electric field evaluator in 1d.
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


#define SOLID1 0.0225
#define SOLID2 0.0875


bool f1( double x, double y, double z )
{
    return( x < SOLID1 );
}


bool f2( double x, double y, double z )
{
    return( x > SOLID2 );
}


void test( void )
{
    Geometry g( MODE_1D, Int3D(11,1,1), Vec3D(0,0,0), 0.01 );
    Solid *s1 = new FuncSolid( f1 );
    g.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( f2 );
    g.set_solid( 8, s2 );
    g.set_boundary( 1, Bound(BOUND_NEUMANN, 0.0) );
    g.set_boundary( 2, Bound(BOUND_DIRICHLET, 0.0) );
    g.set_boundary( 7, Bound(BOUND_DIRICHLET, 0.0) );
    g.set_boundary( 8, Bound(BOUND_DIRICHLET, 0.0) );
    g.build_mesh();
    //g.debug_print();

    EpotProblem p;
    p.construct( g );

    ScalarField epot( g );
    ScalarField scharge( g );
    for( int a = 0; a < g.size(0); a++ )
	scharge(a) = 1.0e-5;

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    EpotEfield ef( g, epot );

    ofstream ostr( "efield.dat" );
    ostr << "# "
	 << setw(10) << "x (m)" << " " 
	 << setw(12) << "potential (V)" << " "
	 << setw(12) << "efield (V/m)" << "\n";
    for( int a = -10; a < 111; a++ ) {
	Vec3D x( g.h()*(g.size(0)-1)*a/99.0, 0.0, 0.0 );
	ostr << setw(12) << x[0] << " " 
	     << setw(12) << epot(x) << " "
	     << setw(12) << ef(x) << "\n";
    }
    ostr << "\n\n";
    ostr.close();

    ostr.open( "efield_solid.dat" );
    ostr << setw(12) << SOLID1 << " "
	 << setw(12) << -10000 << "\n"
	 << setw(12) << SOLID1 << " "
	 << setw(12) << 0 << "\n"
	 << setw(12) << SOLID1 << " "
	 << setw(12) << 10000 << "\n\n\n"
	 << setw(12) << SOLID2 << " "
	 << setw(12) << -10000 << "\n"
	 << setw(12) << SOLID2 << " "
	 << setw(12) << 0 << "\n"
	 << setw(12) << SOLID2 << " "
	 << setw(12) << 10000 << "\n\n\n";
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

