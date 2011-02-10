/*! \file geometry.cpp 
 *  \brief Test geometry definition.
 *
 *  \test Test geometry definition.
 */


#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "geometry.hpp"
#include "func_solid.hpp"
#include "error.hpp"


using namespace std;


bool s1( double x, double y, double z )
{
    return( (x-0.04)*(x-0.04)/(0.032*0.032) + y*y/(0.015*0.015) < 1 );
}


void test1( void )
{
    Geometry g( MODE_2D, Int3D(5,5,1), Vec3D(0,0,0), 0.01 );
    Solid *s = new FuncSolid( s1 );
    g.set_solid( 7, s );
    g.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 2, Bound(BOUND_NEUMANN,   10.0) );
    g.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    g.set_boundary( 4, Bound(BOUND_DIRICHLET, 20.0) );
    g.set_boundary( 7, Bound(BOUND_DIRICHLET, 10.0) );
    g.build_mesh();
    g.debug_print( cout );

    ofstream fstr( "geometry.dat" );
    g.save( fstr );
    fstr.close();
}


void test2( void )
{
    ifstream fstr( "geometry.dat" );
    Geometry g( fstr );
    //g.debug_print();
    int stat = 0;
    try {
	g.inside( Vec3D(0.05,0.05,0.05) );
    } catch( Error ) {
	stat = 1;
    }
    if( stat == 0 )
	exit( 1 );

    /*
    for( int j = 0; j < g.size(1); j++ ) {
	for( int i = 0; i < g.size(0); i++ )
	    cout << setw(3) << (int)g.mesh( i, j );
	cout << "\n";
    }
    */
}


void test( int argc, char **argv )
{
    test1();
    //test2();
}



