/*! \file vecfield_filedef.cpp 
 *  \brief Test VectorField definition from ascii file
 *
 *  \test Test VectorField definition from ascii file
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "vectorfield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"


using namespace std;


void test( int argc, char **argv )
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

