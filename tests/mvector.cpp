/*! \file mvector.cpp 
 * \brief Test non-matrix features of MVector.
 *
 * \test Test non-matrix features of MVector.
 *
 */


#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h>
#define SPM_RANGE_CHECK 1
#include "mvector.hpp"
#include "ibsimutest.hpp"


using namespace std;


void test( int argc, char **argv )
{
    /* Constructors */
    Vector A; 
    if( A.size() != 0 || A.get_data() != NULL )
	throw( ErrorTest( ERROR_LOCATION ) );
    Vector B(5);
    if( B.size() != 5 || B.get_data() == NULL )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( B[0] != 0 || B[1] != 0 || B[2] != 0 || B[3] != 0 || B[4] != 0 )
	throw( ErrorTest( ERROR_LOCATION ) );
    double C_data[5] = {1, -2, 3, -2, 4};
    Vector C(5,C_data); 
    if( C[0] != 1 || C[1] != -2 || C[2] != 3 || C[3] != -2 || C[4] != 4 ) 
	throw( ErrorTest( ERROR_LOCATION ) );
    Vector D(5,2.0); 
    if( D[0] != 2 || D[1] != 2 || D[2] != 2 || D[3] != 2 || D[4] != 2 ) 
	throw( ErrorTest( ERROR_LOCATION ) );
    Vector E = C; 
    if( E[0] != 1 || E[1] != -2 || E[2] != 3 || E[3] != -2 || E[4] != 4 ) 
	throw( ErrorTest( ERROR_LOCATION ) );
    Vector F = C-2.0*D;
    if( F[0] != -3 || F[1] != -6 || F[2] != -1 || F[3] != -6 || F[4] != 0 ) 
	throw( ErrorTest( ERROR_LOCATION ) );


    /* Resizing */
    B.resize( 6 );
    if( B.size() != 6 )
	throw( ErrorTest( ERROR_LOCATION ) );
    B = 3;
    A = B + B;
    if( A.size() != 6 || A[0] != 6 )
	throw( ErrorTest( ERROR_LOCATION ) );
    A.resize( 5 );
    B.resize( 5 );


    /* operator= */
    A = 2;
    if( A[0] != 2 || A[1] != 2 || A[2] != 2 || A[3] != 2 || A[4] != 2 ) 
	throw( ErrorTest( ERROR_LOCATION ) );
    C = A;
    if( C[0] != 2 || C[1] != 2 || C[2] != 2 || C[3] != 2 || C[4] != 2 ) 
	throw( ErrorTest( ERROR_LOCATION ) );


    /* Range checking */
    int stat = 0;
    try {
	A[5] = 2;
    } catch( ErrorRange ) {
	stat = 1;
    }
    if( stat == 0 )
	throw( ErrorTest( ERROR_LOCATION ) );


    /* Dimension checking */
    B.resize( 6 );
    stat = 0;
    try {
	C = A + B;
    } catch( ErrorDim ) {
	stat = 1;
    }
    if( stat == 0 )
	throw( ErrorTest( ERROR_LOCATION ) );
    B.resize( 5 );


    /* Comparison */
    A = 2;
    B = 2*A;
    if( A == B )
	throw( ErrorTest( ERROR_LOCATION ) );
    B = A;
    if( A != B )
	throw( ErrorTest( ERROR_LOCATION ) );


    /* Linear algebra */
    B = 2;
    C = 3;
    D = 0;
    C[2] = -1;
    B[4] = 0;
    D[4] = 1;
    A = 2*B - 1.5*C - D;
    if( A[0] != -0.5 || A[1] != -0.5 || A[2] != 5.5 || A[3] != -0.5 || A[4] != -5.5 )
	throw( ErrorTest( ERROR_LOCATION ) );
    Vector G = B*2 - C*1.5 - D;
    if( G[0] != -0.5 || G[1] != -0.5 || G[2] != 5.5 || G[3] != -0.5 || G[4] != -5.5 )
	throw( ErrorTest( ERROR_LOCATION ) );


    /* Misc vector operations */
    if( norm1(A) != 12.5 )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( fabs( norm2(A) - sqrt(61.25) ) > 0.01 )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( ssqr(A) != 61.25 )
	throw( ErrorTest( ERROR_LOCATION ) );
    A[2] = -1.5;
    A[3] = 1.5;
    if( min(A) != -5.5 )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( min_abs(A) != 0.5 )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( max(A) != 1.5 )
	throw( ErrorTest( ERROR_LOCATION ) );
    if( max_abs(A) != 5.5 )
	throw( ErrorTest( ERROR_LOCATION ) );
    A = 1;
    B = 2;
    if( dot_prod( A, B ) != 10.0 )
	throw( ErrorTest( ERROR_LOCATION ) );


    /* Data operations */
    A = 1;
    B = 2;
    swap( A, B );
    if( A[0] != 2 || B[0] != 1 )
	throw( ErrorTest( ERROR_LOCATION ) );
    A.merge( B );
    if( A[0] != 1 || B.size() != 0 )
	throw( ErrorTest( ERROR_LOCATION ) );


    /* Initialization */
    A = 2;
    if( A[0] != 2 || A[1] != 2 || A[2] != 2 || A[3] != 2 || A[4] != 2 )
	throw( ErrorTest( ERROR_LOCATION ) );
}



