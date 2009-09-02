/*! \file efield_extrpl.cpp 
 *  \brief Test for electric field evaluator extrapolation in 1d.
 *
 *  \test Test produces a parabel shape potential between dirichlet
 *  boundaries at x=0 and x=1. Electric field is evaluated from x=-1.5
 *  to x=2.5 with different boundary extrapolation settings. No
 *  automatic verification of test results.
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


void test( void )
{
    Geometry g( MODE_1D, Int3D(11,1,1), Vec3D(0,0,0), 0.1 );
    g.set_boundary( 1, Bound(BOUND_DIRICHLET, 0.0) );
    g.set_boundary( 2, Bound(BOUND_DIRICHLET, 0.0) );
    g.build_mesh();
    //g.debug_print();

    EpotProblem p;
    p.construct( g );

    ScalarField epot( g );
    ScalarField scharge( g );
    for( int a = 0; a < g.size(0); a++ )
	scharge(a) = 1.0e-8;

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    EpotEfield ef( g, epot );
    efield_extrpl_e extrpl[6] = { EFIELD_MIRROR, EFIELD_MIRROR,
				  EFIELD_MIRROR, EFIELD_MIRROR,
				  EFIELD_MIRROR, EFIELD_MIRROR };
    ef.set_extrapolation( extrpl );

    ofstream ostr( "efield_extrpl.dat" );
    ostr << "# "
	 << setw(10) << "x (m)" << " " 
	 << setw(12) << "potential (V)" << " "
	 << setw(12) << "efield (V/m)" << "\n";
    for( int a = -30; a < 50; a++ ) {
	Vec3D x( g.h()*(g.size(0)-1)*a/20.0, 0.0, 0.0 );
	ostr << setw(12) << x[0] << " " 
	     << setw(12) << epot(x) << " "
	     << setw(12) << ef(x)[0] << "\n";
    }
    ostr << "\n\n";


    efield_extrpl_e extrpl2[6] = { EFIELD_ZERO, EFIELD_ZERO,
				   EFIELD_ZERO, EFIELD_ZERO,
				   EFIELD_ZERO, EFIELD_ZERO };
    ef.set_extrapolation( extrpl2 );

    for( int a = -30; a < 50; a++ ) {
	Vec3D x( g.h()*(g.size(0)-1)*a/20.0, 0.0, 0.0 );
	ostr << setw(12) << x[0] << " " 
	     << setw(12) << epot(x) << " "
	     << setw(12) << ef(x)[0] << "\n";
    }
    ostr << "\n\n";


    efield_extrpl_e extrpl3[6] = { EFIELD_EXTRAPOLATE, EFIELD_EXTRAPOLATE,
				   EFIELD_EXTRAPOLATE, EFIELD_EXTRAPOLATE,
				   EFIELD_EXTRAPOLATE, EFIELD_EXTRAPOLATE };
    ef.set_extrapolation( extrpl3 );

    for( int a = -30; a < 50; a++ ) {
	Vec3D x( g.h()*(g.size(0)-1)*a/20.0, 0.0, 0.0 );
	ostr << setw(12) << x[0] << " " 
	     << setw(12) << epot(x) << " "
	     << setw(12) << ef(x)[0] << "\n";
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

