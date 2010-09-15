/*! \file vecfield_filedef.cpp 
 *  \test Test VectorField definition from ascii file
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "vectorfield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"


using namespace std;


void test( void )
{
    bool fout[3] = {false, false, true};
    char *srcdir = getenv( "srcdir" );
    char *srcdir_default = ".";
    if( srcdir == NULL )
	srcdir = srcdir_default;
    string infile = to_string(srcdir) + "/vecfield_filedef.in";

    VectorField f( MODE_2D, fout, 1.0, 1.0, infile );
    //f.debug_print();
}


int main( void )
{
    try {
	ibsimu.set_verbose_output( 0 );
	test();
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}



