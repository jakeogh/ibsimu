/*! \file efieldtest_2d.cpp 
 *  \brief Test electric field interpolation in 2d.
 *
 *  \test Test electric field interpolation in 2d.
 *
 */


#include <fstream>
#include <iomanip>
#include "epot_efield.hpp"
#include "error.hpp"


using namespace std;


void test( void )
{
    Geometry g( MODE_2D, Int3D(3,3,1), Vec3D(0,0,0), 1.0 );
    ScalarField epot( g );
    epot(0,0) = 1.0;
    epot(1,0) = 2.0;
    epot(2,0) = 1.0;
    epot(0,1) = 1.0;
    epot(1,1) = 2.0;
    epot(2,1) = 1.0;
    epot(0,2) = 1.0;
    epot(1,2) = 2.0;
    epot(2,2) = 1.0;
    EpotEfield ef( g, epot );

    ofstream ostr( "efieldtest_2d.dat" );
    for( int a = 0; a < 5; a++ ) {
	for( int b = 0; b < 5; b++ ) {
	    Vec3D x( g.h()*(g.size(0)-1)*a/4.0, 
		     g.h()*(g.size(0)-1)*b/4.0, 
		     0.0 );
	    ostr << setw(12) << x << " " 
		 << setw(12) << epot(x) << " "
		 << setw(12) << ef(x) << "\n";
	}
    }
    ostr << "\n\n";
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

