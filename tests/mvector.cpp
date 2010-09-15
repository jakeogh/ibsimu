/* MVector test -- Test non-matrix features of MVector.
 *
 *
 */


#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h>
#define SPM_RANGE_CHECK 1
#include "mvector.hpp"


using namespace std;


#define ERROR()							    \
    {								    \
	cout << "Error at " << __FILE__ << ":" << __LINE__ << "\n"; \
	return( 1 );					    \
    }


int main( void )
{
    /* Constructors */
    Vector A; 
    if( A.size() != 0 || A.get_data() != NULL )
	ERROR();
    Vector B(5);
    if( B.size() != 5 || B.get_data() == NULL )
	ERROR();
    if( B[0] != 0 || B[1] != 0 | B[2] != 0 || B[3] != 0 || B[4] != 0 )
	ERROR();
    double C_data[5] = {1, -2, 3, -2, 4};
    Vector C(5,C_data); 
    if( C[0] != 1 || C[1] != -2 || C[2] != 3 || C[3] != -2 || C[4] != 4 ) 
	ERROR();
    Vector D(5,2.0); 
    if( D[0] != 2 || D[1] != 2 || D[2] != 2 || D[3] != 2 || D[4] != 2 ) 
	ERROR();
    Vector E = C; 
    if( E[0] != 1 || E[1] != -2 || E[2] != 3 || E[3] != -2 || E[4] != 4 ) 
	ERROR();
    Vector F = C-2.0*D;
    if( F[0] != -3 || F[1] != -6 || F[2] != -1 || F[3] != -6 || F[4] != 0 ) 
	ERROR();


    /* Resizing */
    B.resize( 6 );
    if( B.size() != 6 )
	ERROR();
    B = 3;
    A = B + B;
    if( A.size() != 6 || A[0] != 6 )
	ERROR();
    A.resize( 5 );
    B.resize( 5 );


    /* operator= */
    A = 2;
    if( A[0] != 2 || A[1] != 2 || A[2] != 2 || A[3] != 2 || A[4] != 2 ) 
	ERROR();
    C = A;
    if( C[0] != 2 || C[1] != 2 || C[2] != 2 || C[3] != 2 || C[4] != 2 ) 
	ERROR();


    /* Range checking */
    int stat = 0;
    try {
	A[5] = 2;
    } catch( ErrorRange ) {
	stat = 1;
    }
    if( stat == 0 )
	ERROR();


    /* Dimension checking */
    B.resize( 6 );
    stat = 0;
    try {
	C = A + B;
    } catch( ErrorDim ) {
	stat = 1;
    }
    if( stat == 0 )
	ERROR();
    B.resize( 5 );


    /* Comparison */
    A = 2;
    B = 2*A;
    if( A == B )
	ERROR();
    B = A;
    if( A != B )
	ERROR();


    /* Linear algebra */
    B = 2;
    C = 3;
    D = 0;
    C[2] = -1;
    B[4] = 0;
    D[4] = 1;
    A = 2*B - 1.5*C - D;
    if( A[0] != -0.5 || A[1] != -0.5 || A[2] != 5.5 || A[3] != -0.5 || A[4] != -5.5 )
	ERROR();
    Vector G = B*2 - C*1.5 - D;
    if( G[0] != -0.5 || G[1] != -0.5 || G[2] != 5.5 || G[3] != -0.5 || G[4] != -5.5 )
	ERROR();


    /* Misc vector operations */
    if( norm1(A) != 12.5 )
	ERROR();
    if( fabs( norm2(A) - sqrt(61.25) ) > 0.01 )
	ERROR();
    if( ssqr(A) != 61.25 )
	ERROR();
    A[2] = -1.5;
    A[3] = 1.5;
    if( min(A) != -5.5 )
	ERROR();
    if( min_abs(A) != 0.5 )
	ERROR();
    if( max(A) != 1.5 )
	ERROR();
    if( max_abs(A) != 5.5 )
	ERROR();
    A = 1;
    B = 2;
    if( dot_prod( A, B ) != 10.0 )
	ERROR();


    /* Data operations */
    A = 1;
    B = 2;
    swap( A, B );
    if( A[0] != 2 || B[0] != 1 )
	ERROR();
    A.merge( B );
    if( A[0] != 1 || B.size() != 0 )
	ERROR();


    /* Initialization */
    A = 2;
    if( A[0] != 2 || A[1] != 2 || A[2] != 2 || A[3] != 2 || A[4] != 2 )
	ERROR();

    
    return( 0 );
}



