/*! \file vecfield_filedef.cpp 
 *  \test Test VectorField definition from ascii file
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "vectorfield.hpp"
#include "verbose.hpp"
#include "error.hpp"


using namespace std;


void test( void )
{
    bool fout[3] = {false, false, true};
    VectorField f( MODE_2D, fout, 1.0, 1.0, "vecfield_filedef.in" );
    //f.debug_print();
}


int main( void )
{
    try {
	verbose_output = 0;
	test();
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

