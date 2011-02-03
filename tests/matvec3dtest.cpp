/*! \file matvec3dtest.cpp 
 *  \brief Basic tests for 3D vectors and matrices
 *
 *  \test Basic tests for 3D vectors and matrices
 *
 */


#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "vec3d.hpp"
#include "mat3d.hpp"


using namespace std;


//#define PRINT 1

#define ERROR()							    \
    {								    \
	cout << "Error at " << __FILE__ << ":" << __LINE__ << "\n"; \
	exit( 1 );						    \
    }


int test( int argc, char **argv )
{
    /* Vector ssqr, norm2, normalize
     */
    Vec3D a(2,-1,2);
    if( fabs(a.ssqr()-9) > 1e-6 )
	ERROR();
    if( fabs(a.norm2()-3) > 1e-6 )
	ERROR();
    if( fabs(norm2(a)-3) > 1e-6 )
	ERROR();
    a.normalize();
    if( norm2(a-Vec3D(2.0/3.0,-1.0/3.0,2.0/3.0)) > 1e-6 )
	ERROR();

    /* Matrix determinant, matrix inversion,
     * vector-matrix multiplication
     */
    Vec3D x( 1, 2, 3 );
    Mat3D m(  1, -2,  3, 
	      1,  2,  2,
	     -1,  4, -2 );
    double mdet = m.determinant();
#ifdef PRINT
    cout << "x = \n" << x << "\n";
    cout << "m = \n" << m << "\n";
    cout << "mdet = " << mdet << "\n";
#endif
    if( fabs(mdet-6.0) > 1e-6 )
	ERROR();
    Vec3D b = m*x;
#ifdef PRINT
    cout << "b = m*x = \n" << b << "\n";
#endif
    if( norm2(b-Vec3D(6,11,1)) > 1e-6 )
	ERROR();
    Mat3D minv = m.inverse();
    Vec3D xnew = minv*b;
#ifdef PRINT
    cout << "minv = \n" << minv << "\n";
    cout << "xnew = minv*b = \n" << xnew << "\n";
#endif
    if( norm2(xnew-x) > 1e-6 )
	ERROR();

    /* Dot product
     */
    double dot = x*b;
#ifdef PRINT
    cout << "x*b = " << dot << "\n";
#endif
    if( fabs(dot-31.0) > 1e-6 )
	ERROR();

    /* Cross product
     */
    Vec3D y = cross(x,b);
#ifdef PRINT
    cout << "y = cross(x,b) = \n" << y << "\n";
#endif
    if( norm2(y-Vec3D(-31,17,-1)) > 1e-6 )
	ERROR();

    return( 0 );
}
