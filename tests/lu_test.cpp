/* LU test -- Test ILU matrix decomposition and solvers
 *
 *
 */


#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h>
#define SPM_RANGE_CHECK 1
#include "crowmatrix.hpp"
#include "ccolmatrix.hpp"
#include "mvector.hpp"
#include "ilu0_precond.hpp"


using namespace std;


#define ERROR()							    \
    {								    \
	cout << "Error at " << __FILE__ << ":" << __LINE__ << "\n"; \
	return( 1 );					    \
    }


void init_fdm_problem( Matrix &A, Vector &B, double v1, double v2 )
{
    A.resize( 15, 15 );
    B.resize( 15 );

    A.set(0,0) = -1;
    A.set(0,3) = 1;
    B(0) = 0.0;

    A.set(1,1) = -1;
    A.set(1,4) = 1;
    B(1) = 0.0;

    A.set(2,2) = -1;
    A.set(2,5) = 1;
    B(2) = 0.0;

    A.set(3,0) = 1;
    A.set(3,3) = -4;
    A.set(3,4) = 1;
    A.set(3,6) = 1;
    B(3) = -v1;

    A.set(4,1) = 1;
    A.set(4,3) = 1;
    A.set(4,4) = -4;
    A.set(4,5) = 1;
    A.set(4,7) = 1;
    B(4) = 0.0;

    A.set(5,2) = 1;
    A.set(5,4) = 1;
    A.set(5,5) = -4;
    A.set(5,8) = 1;
    B(5) = -v2;

    A.set(6,3) = 1;
    A.set(6,6) = -4;
    A.set(6,7) = 1;
    A.set(6,9) = 1;
    B(6) = -v1;

    A.set(7,4) = 1;
    A.set(7,6) = 1;
    A.set(7,7) = -4;
    A.set(7,8) = 1;
    A.set(7,10) = 1;
    B(7) = 0.0;

    A.set(8,5) = 1;
    A.set(8,7) = 1;
    A.set(8,8) = -4;
    A.set(8,11) = 1;
    B(8) = -v2;

    A.set(9,6) = 1;
    A.set(9,9) = -4;
    A.set(9,10) = 1;
    A.set(9,12) = 1;
    B(9) = -v1;

    A.set(10,7) = 1;
    A.set(10,9) = 1;
    A.set(10,10) = -4;
    A.set(10,11) = 1;
    A.set(10,13) = 1;
    B(10) = -1.129e3;  // 1e-4 C/m^3 on a 1 cm mesh

    A.set(11,8) = 1;
    A.set(11,10) = 1;
    A.set(11,11) = -4;
    A.set(11,14) = 1;
    B(11) = -v2;

    A.set(12,9) = 1;
    A.set(12,12) = -1;
    B(12) = 0.0;

    A.set(13,10) = 1;
    A.set(13,13) = -1;
    B(13) = 0.0;

    A.set(14,11) = 1;
    A.set(14,14) = -1;
    B(14) = 0.0;
}


int main( void )
{
    CColMatrix Ac;
    CRowMatrix Ar;
    Vector B, X, Y;

    double X_check_data[15] = { 54.1396,
				66.1492,
				10.3966,
				54.1396,
				66.1492,
				10.3966,
				86.2697,
				162.668,
				36.0411,
				140.321,
				504.379,
				131.559,
				5,
				466.062,
				131.559 };
    Vector X_check( 15, X_check_data );
    double Y_check_data[15] = { 0,
				0,
				0,
				-10,
				-3.33333,
				8.75,
				-13.3333,
				-4.88636,
				11.8757,
				-13.6364,
				-1134.12,
				-316.742,
				-3.65854,
				-330.318,
				-92.9515 };
    Vector Y_check( 15, Y_check_data );

    // CColMatrix 
    init_fdm_problem( Ac, B, 10.0, -10.0 );
    ILU0_Precond pc = ILU0_Precond( Ac );
    const Matrix *Lc = pc.get_L();
    const Matrix *Uc = pc.get_U();
    Lc->lower_unit_solve( Y, B );
    Uc->upper_diag_solve( X, Y );

    if( max_abs(Y-Y_check) > 0.1 )
	ERROR();
    if( max_abs(X-X_check) > 0.1 )
	ERROR();

    // CRowMatrix
    init_fdm_problem( Ar, B, 10.0, -10.0 );
    ILU0_Precond pr = ILU0_Precond( Ar );
    const Matrix *Lr = pr.get_L();
    const Matrix *Ur = pr.get_U();
    Lr->lower_unit_solve( Y, B );
    Ur->upper_diag_solve( X, Y );

    if( max_abs(Y-Y_check) > 0.1 )
	ERROR();
    if( max_abs(X-X_check) > 0.1 )
	ERROR();

    return( 0 );
}

