/*! \file vecfield_filedef.cpp 
 *  \brief Test VectorField definition from ascii file.
 *
 *  \test Test VectorField definition from ascii file.
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "meshvectorfield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"


using namespace std;


void test( int argc, char **argv )
{
    const char *srcdir = getenv( "srcdir" );
    const char *srcdir_default = ".";
    if( srcdir == NULL )
	srcdir = srcdir_default;
    string infile = to_string(srcdir) + "/vecfield_filedef.in";

    bool fout[3] = {false, false, true};
    MeshVectorField f( MODE_2D, fout, 1.0e-3, 1000.0, infile );
    f.debug_print( cout );
}

