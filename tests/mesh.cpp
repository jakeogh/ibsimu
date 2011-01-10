/*! \file mesh.cpp 
 * \brief Test mesh class.
 *
 * \test Test mesh class.
 *
 */



#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "geometry.hpp"
#include "error.hpp"


using namespace std;


void test( int argc, char **argv )
{
    Mesh m1( MODE_2D, Int3D(3,4,5), Vec3D(0.0,-1.0,2.0), 1.0 );
    Mesh m2( m1 );
    if( !(m1 == m2) ) {
	std::cout << "Error\n";
	exit( 1 );
    }
}


int main( int argc, char **argv )
{
    try {
	test( argc, argv );
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	e.print_trace( cout );
	exit( 1 );
    }

    return( 0 );
}



