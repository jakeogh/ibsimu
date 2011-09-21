/*! \file vecfield_rotate.cpp 
 *  \brief Test VectorField rotate functions
 *
 *  \test Test VectorField rotate functions
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "meshvectorfield.hpp"
#include "ibsimutest.hpp"
#include "error.hpp"


using namespace std;


void test( int argc, char **argv )
{
    bool fout[3] = {true, true, false};
    Int3D size(2,3,1);
    Vec3D origo(2.2,1.1,0);
    MeshVectorField f( MODE_2D, fout, size, origo, 1.0 );
    f.set( 0, 0, Vec3D(  0.1,  0.4, 0 ) );
    f.set( 1, 0, Vec3D(  0.2,  1.5, 0 ) );
    f.set( 0, 1, Vec3D(  0.3,  1.0, 0 ) );
    f.set( 1, 1, Vec3D( -0.4, -1.5, 0 ) );
    f.set( 0, 2, Vec3D(  0.5,  4.0, 0 ) );
    f.set( 1, 2, Vec3D(  0.6,  3.0, 0 ) );

    if( norm2( Vec3D(0.1,0.4,0)-f(Vec3D(2.2,1.1,0)) ) > 0.01 )
	throw( ErrorTest( ERROR_LOCATION ) );

    f.rotate_z( 90 );

    if( norm2( Vec3D(-0.4,0.1,0)-f(Vec3D(-1.1,2.2,0)) ) > 0.01 )
	throw( ErrorTest( ERROR_LOCATION ) );

    f.rotate_z( -90 );

    if( norm2( Vec3D(0.1,0.4,0)-f(Vec3D(2.2,1.1,0)) ) > 0.01 )
	throw( ErrorTest( ERROR_LOCATION ) );

    f.rotate_z( 180 );

    if( norm2( Vec3D(-0.1,-0.4,0)-f(Vec3D(-2.2,-1.1,0)) ) > 0.01 )
	throw( ErrorTest( ERROR_LOCATION ) );

    f.scale( 0.5 );

    if( norm2( Vec3D(-0.1,-0.4,0)-f(Vec3D(-1.1,-0.55,0)) ) > 0.01 )
	throw( ErrorTest( ERROR_LOCATION ) );

    f.scale( -2.0 );

    if( norm2( Vec3D(0.1,0.4,0)-f(Vec3D(2.2,1.1,0)) ) > 0.01 )
	throw( ErrorTest( ERROR_LOCATION ) );
}




