/*! \file matrixsolver.cpp 
 *  \test Test internal matrix solvers.
 */


#include <cstdlib>
#include "bicgstab.hpp"
#include "crowmatrix.hpp"
#include "ccolmatrix.hpp"
#include "mvector.hpp"
#include "empty_precond.hpp"
#include "diag_precond.hpp"
#include "ilu0_precond.hpp"
#include "gs_solver.hpp"
#include "hbio.hpp"


using namespace std;


void test( void )
{
    Vector rhs, sol[4];
    CColMatrix cmat;
    CRowMatrix mat;
    HBIO hb;
    double x, w;
    int imax[4];
    double eps[4];
    bool failed = false;

    char *srcdir = getenv( "srcdir" );
    char *srcdir_default = ".";
    if( srcdir == NULL )
	srcdir = srcdir_default;
    string filename = to_string(srcdir) + "\matrixsolver.hb";

    hb.read( filename );
    hb.get_matrix( cmat );
    hb.get_rhs_vector( rhs );
    mat = cmat;

    for( int a = 0; a < 4; a++ ) {
	sol[a].resize( rhs.size() );
	sol[a] = 0;
    }

    //std::cout << "Solving BiCGSTAB with empty preconditioner\n";
    eps[0] = 1.0e-9;
    imax[0] = 100000;
    if( !bicgstab( cmat, rhs, sol[0], Empty_Precond(), imax[0], eps[0] ) ) {
	cout << "BiCGSTAB with empty preconditioner failed.\n";
	exit( 1 );
    }

    //std::cout << "Solving BiCGSTAB with diagonal preconditioner\n";
    eps[1] = 1.0e-9;
    imax[1] = 100000;
    if( !bicgstab( cmat, rhs, sol[1], Diag_Precond(cmat), imax[1], eps[1] ) ) {
	cout << "BiCGSTAB with diagonal preconditioner failed.\n";
	exit( 1 );
    }

    //std::cout << "Solving BiCGSTAB with ILU0 preconditioner\n";
    eps[2] = 1.0e-9;
    imax[2] = 100000;
    if( !bicgstab( cmat, rhs, sol[2], ILU0_Precond(cmat), imax[2], eps[2] ) ) {
	cout << "BiCGSTAB with ILU0 preconditioner failed.\n";
	exit( 1 );
    }

    //std::cout << "Solving Gauss-Seidel\n";
    w = 1.66;
    eps[3] = 1.0e-9;
    imax[3] = 100000;
    if( !GSSolver::gauss_seidel( mat, rhs, sol[3], imax[3], eps[3], w ) ) {
	cout << "Gauss-Seidel solver failed.\n";
	exit( 1 );
    }

    /* Compare results */
    for( int a = 1; a < 4; a++ ) {
	x = max_abs(sol[0]-sol[a]);
	if( x >  1.0e-4 )
	    failed = true;
    }

    if( failed ) {
	/* Print diagnostics if failed */
	for( int a = 0; a < 4; a++ ) {
	    std::cout << "Solver " << a 
		      << ": imax = " << imax[a] 
		      << " eps = " << eps[a] 
		      << " diff = " << max_abs(sol[0]-sol[a])
		      << "\n";
	}

	exit( 1 );
    }
}


int main( void )
{
    try {
	test();
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

