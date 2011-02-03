#include <iostream>
#include "ibsimutest.hpp"
#include "ibsimu.hpp"
#include "error.hpp"


/* Main function for tests
 */
int main( int argc, char **argv )
{
    try {
	ibsimu.set_verbose_output( 0 );
	ibsimu.set_thread_count( 4 );
	test( argc, argv );
    } catch( Error e ) {
	e.print_error_message( std::cerr );
	exit( 1 );
    }

    return( 0 );
}
