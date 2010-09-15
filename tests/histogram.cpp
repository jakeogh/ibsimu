/*! \file histogram.cpp 
 *  \test Test histogram.
 *
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include "histogram.hpp"
#include "error.hpp"


void test( int *argc, char ***argv )
{
    double range[4] = {0, 0, 1, 1};
    Histogram2D histo( 4, 4, range );

    histo.accumulate_linear( 1.0, 1.0/3.0, 1.0 );
    histo.accumulate_linear( 0.5, 0.5, 1.0 );
    histo.accumulate_linear( 3.0, 3.0, 1.0 );

    if( fabs(histo(0,0) - 0.0) > 1.0e-6 ||
	fabs(histo(0,1) - 0.0) > 1.0e-6 ||
	fabs(histo(0,2) - 0.0) > 1.0e-6 ||
	fabs(histo(0,3) - 0.0) > 1.0e-6 ||

	fabs(histo(1,0) - 0.0) > 1.0e-6 ||
	fabs(histo(1,1) - 0.25) > 1.0e-6 ||
	fabs(histo(1,2) - 0.25) > 1.0e-6 ||
	fabs(histo(1,3) - 0.0) > 1.0e-6 ||

	fabs(histo(2,0) - 0.0) > 1.0e-6 ||
	fabs(histo(2,1) - 0.25) > 1.0e-6 ||
	fabs(histo(2,2) - 0.25) > 1.0e-6 ||
	fabs(histo(2,3) - 0.0) > 1.0e-6 ||

	fabs(histo(3,0) - 0.0) > 1.0e-6 ||
	fabs(histo(3,1) - 1.0) > 1.0e-6 ||
	fabs(histo(3,2) - 0.0) > 1.0e-6 ||
	fabs(histo(3,3) - 0.0) > 1.0e-6 ) 
    {
	std::cout << "Histogram value not expected\n\n";
	for( size_t j = 0; j < 4; j++ ) {
	    for( size_t i = 0; i < 4; i++ ) {
		std::cout << std::setw(12) << histo(i,j) << " ";
	    }
	    std::cout << "\n";
	}

	exit( 1 );
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



