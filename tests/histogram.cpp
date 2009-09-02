/*! \file histogram.cpp 
 *  \brief Test histogram.
 *
 *  \test Test histogram.
 *
 */

#include <iostream>
#include <iomanip>
#include <vector>

#include "histogram.hpp"
#include "error.hpp"


void test( int *argc, char ***argv )
{
    double range[4] = {0, 0, 1, 1};
    Histogram2D histo( 4, 4, range );

    histo.accumulate_linear( 0.5, 0.5, 1.0 );

    for( int j = histo.m()-1; j >= 0; j-- ) {
	for( int i = 0; i < (int)histo.n(); i++ ) {
	    std::cout << histo( i, j ) << " ";
	}
	std::cout << "\n";
    }
}


int main( int argc, char **argv )
{
    try {
	test( &argc, &argv );
    } catch ( Error e ) {
	std::cout << "Error in " << e._loc._file << ":" << e._loc._line 
		  << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

