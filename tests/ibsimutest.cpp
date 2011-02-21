/*! \file ibsimutest.cpp
 *  \brief Common functions for tests.
 */


#include <iostream>
#include "ibsimutest.hpp"
#include "ibsimu.hpp"


ErrorTest::ErrorTest( const ErrorLocation &loc, const std::string &str )
    : Error(loc)
{
    _error_str = "Test failure: " + str;
}


/* Main function for tests
 */
int main( int argc, char **argv )
{
    try {
	ibsimu.set_verbose_output( 1 );
	ibsimu.set_thread_count( 2 );
	test( argc, argv );
    } catch( Error e ) {
	e.print_error_message( std::cerr );
	exit( 1 );
    }

    return( 0 );
}
