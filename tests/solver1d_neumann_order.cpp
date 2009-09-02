/*! \file solver1d_neumann_order.cpp 
 *  \brief Test solver with a 1d problem with different Neumann boundary orders.
 *
 *  \test  Test solver with a 1d problem with different Neumann boundary orders.
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


void test( void )
{
    Geometry g( MODE_1D, Int3D(11,1,1), Vec3D(0,0,0), 0.01 );
    g.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 2, Bound(BOUND_NEUMANN,   10.0) );
    g.build_mesh();

    EpotProblem p;
    p.set_neumann_order( 1 );
    p.construct( g );

    ScalarField epot( g );
    ScalarField scharge( g );
    for( int a = 0; a < scharge.size(0); a++ )
	scharge(a) = 1.0e-4;

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    ofstream ostr;
    ostr.open( "solver1d_neumann_order1.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "potential (V)" << "\n";
    for( int a = 0; a < g.size(0); a++ ) {
	ostr << setw(14) << a*g.h() << " " 
	     << setw(14) << epot(a) << "\n";
    }
    ostr << "\n\n";
    ostr.close();

    // *********************************

    p.set_neumann_order( 2 );
    p.construct( g );
    p.solve( epot, scharge );

    ostr.open( "solver1d_neumann_order2.dat" );
    ostr << "# "
	 << setw(12) << "x (m)" << " " 
	 << setw(14) << "potential (V)" << "\n";
    for( int a = 0; a < g.size(0); a++ ) {
	ostr << setw(14) << a*g.h() << " " 
	     << setw(14) << epot(a) << "\n";
    }
    ostr << "\n\n";
    ostr.close();
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

