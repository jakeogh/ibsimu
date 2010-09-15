/*! \file sclfield.cpp
 *  \test Test for ScalarField class
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "scalarfield.hpp"
#include "error.hpp"


using namespace std;


void test( void )
{
    ScalarField f;
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



