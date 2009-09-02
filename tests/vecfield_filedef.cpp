/*! \file vecfield_filedef.cpp 
 *  \brief Test VectorField definition from ascii file
 *
 *  \test Test VectorField definition from ascii file
 */


#include <iostream>
#include <iomanip>
#include "vectorfield.hpp"
#include "verbose.hpp"
#include "error.hpp"


using namespace std;


void test( void )
{
    verbose_output = 1;

    bool fout[3] = {false, false, true};
    VectorField f( MODE_2D, fout, 1.0, 1.0, "vecfield_filedef.dat" );
    f.debug_print();
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

