/*! \file solver1d_solid_neumann_conflict.cpp 
 *  \brief Test solver with a 1d problem with Neumann and smooth 
 *  solid conflicting
 *
 *  \test  Test solver with a 1d problem with Neumann and smooth 
 *  solid conflicting
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
#include "verbose.hpp"


using namespace std;


bool solid1( double x, double y, double z )
{
    // Leave nodes at 0.09 and 0.1 free = OK
    //return( x >= 0.055 && x <= 0.085 );

    // Leave only node at 0.1 free
    return( x >= 0.055 && x <= 0.095 );
}


void test( void )
{
    verbose_output = 1;
    Geometry g( MODE_1D, Int3D(11,1,1), Vec3D(0,0,0), 0.01 );
    Solid *s1 = new FuncSolid( solid1 );
    g.set_solid( 7, s1 );
    g.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    g.set_boundary( 2, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 7, Bound(BOUND_DIRICHLET, 10.0) );
    g.build_mesh();
    g.debug_print();

    EpotProblem p;
    p.set_neumann_order( 2 );
    p.enable_smooth_solids( true );
    p.construct( g );
    p.debug_print();

    ScalarField epot( g );
    ScalarField scharge( g );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    ofstream ostr( "solver1d_solid_neumann_conflict.dat" );
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

