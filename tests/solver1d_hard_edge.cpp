/*! \file solver1d_hard_edge.cpp 
 *  \brief Test solver with a 1d problem with hard edges.
 *
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


double x1 = 0.095;


bool s1( double x, double y, double z )
{
    return( x > x1 );
}


void test( void )
{
    Geometry g( MODE_1D, Int3D(11,1,1), Vec3D(0,0,0), 0.01 );
    Solid *s = new FuncSolid( s1 );
    g.set_solid( 7, s );
    g.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 2, Bound(BOUND_DIRICHLET, 10.0) );
    g.set_boundary( 7, Bound(BOUND_DIRICHLET, 10.0) );
    g.build_mesh();

    EpotProblem p;
    p.enable_smooth_solids( false );
    p.construct( g );

    g.debug_print();
    p.debug_print();

    ScalarField epot( g );
    ScalarField scharge( g );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );
    
    ofstream ostr( "solver1d_hard_edge.dat" );
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

