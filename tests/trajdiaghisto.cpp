/*! \file trajdiaghisto.cpp 
 *  \brief Test particle trajectory diagnostic histogram
 *
 *  \test Test particle trajectory diagnostic histogram
 *
 */


#include <fstream>
#include <iomanip>
#include "fileplot.hpp"
#include "geometry.hpp"
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "epot_efield.hpp"
#include "particles.hpp"
#include "error.hpp"
#include "verbose.hpp"


using namespace std;


#define ERROR()							    \
    {								    \
	cout << "Error at " << __FILE__ << ":" << __LINE__ << "\n"; \
	exit( 1 );						    \
    }


void test( void )
{
    TrajectoryDiagnosticColumn c1( DIAG_X );
    TrajectoryDiagnosticColumn c2( DIAG_XP );

    c1.add_data( -5.0 );
    c2.add_data( -5.0 );

    c1.add_data(  5.0 );
    c2.add_data(  5.0 );

    c1.add_data(  0.0 );
    c2.add_data(  0.0 );

    Histogram2D h( 11, 11, c1, c2 );

    for( int j = 10; j >= 0; j-- ) {
	for( int i = 0; i < 11; i++ ) {
	    cout << setw(12) << h(i,j) << " ";
	}
	cout << "\n";
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

