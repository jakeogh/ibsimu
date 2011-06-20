/*! \file mssolver1d_neumann.cpp 
 *  \brief Test multi grid solver with a 1d problem with dirichlet+neumann boundaries
 *
 *  \test Test multi grid solver with a 1d problem with dirichlet+neumann boundaries
 */


#include <fstream>
#include <iomanip>
#include "epot_gssolver.hpp"
#include "epot_mgsolver.hpp"
#include "geometry.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "constants.hpp"


using namespace std;


const double rho = 1.0e-4;
const double A = rho/EPSILON0*0.04;
const double B = 0.0;


double phi( double x )
{
    return( -rho/(2.0*EPSILON0)*x*x + A*x + B );
}


void test( int argc, char **argv )
{
    bool err = false;
    //Geometry geom( MODE_1D, Int3D(3,1,1),   Vec3D(0,0,0), 0.02 );      // levels 1
    //Geometry geom( MODE_1D, Int3D(5,1,1),   Vec3D(0,0,0), 0.01 );      // levels 2
    //Geometry geom( MODE_1D, Int3D(9,1,1),   Vec3D(0,0,0), 0.005 );     // levels 3
    //Geometry geom( MODE_1D, Int3D(17,1,1),  Vec3D(0,0,0), 0.0025 );    // levels 4
    //Geometry geom( MODE_1D, Int3D(33,1,1),  Vec3D(0,0,0), 0.00125 );   // levels 5
    //Geometry geom( MODE_1D, Int3D(65,1,1),  Vec3D(0,0,0), 0.000625 );  // levels 6
    //Geometry geom( MODE_1D, Int3D(129,1,1), Vec3D(0,0,0), 0.0003125 ); // levels 7
    //Geometry geom( MODE_1D, Int3D(257,1,1), Vec3D(0,0,0), 1.5625e-4 ); // levels 8
    Geometry geom( MODE_1D, Int3D(513,1,1), Vec3D(0,0,0), 7.8125e-5 ); // levels 9
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,    0.0) );
    geom.build_mesh();

    MeshScalarField scharge( geom );
    for( size_t i = 0; i < scharge.size(0); i++ ) {
	scharge(i) = rho;
    }

    EpotMGSolver solver( geom );
    MeshScalarField epot( geom );
    solver.set_levels( 8 );
    solver.set_gamma( 1 );
    solver.set_neumann_order( 2 );
    solver.solve( epot, scharge );
    
    if( true ) {
	std::cout << "epot (MG):\n";
	double eps = 0.0;
	size_t loc = 0;
	for( size_t i = 0; i < geom.size(0); i++ ) {
	    std::cout << std::setw(8) << epot(i) << " ";
	    double dif = fabs( epot(i) - phi(i*geom.h()) );
	    if( dif > eps ) {
		eps = dif;
		loc = i;
	    }
	}
	std::cout << "\n";
	std::cout << "Maximum difference to analytic: " << eps << " at " << loc << "\n";
    }

    EpotGSSolver solver2( geom );
    MeshScalarField epot2( geom );
    solver2.set_neumann_order( 2 );
    solver2.set_imax( 100000 );
    solver2.solve( epot2, scharge );

    if( true ) {
	std::cout << "epot2 (GS):\n";
	double eps = 0.0;
	size_t loc = 0;
	for( size_t i = 0; i < geom.size(0); i++ ) {
	    std::cout << std::setw(8) << epot2(i) << " ";
	    double dif = fabs( epot2(i) - phi(i*geom.h()) );
	    if( dif > eps ) {
		eps = dif;
		loc = i;
	    }
	}
	std::cout << "\n";
	std::cout << "Maximum difference to analytic: " << eps << " at " << loc << "\n";
    }

    if( true ) {
	double eps = 0.0;
	size_t loc = 0;
	for( size_t i = 0; i < geom.size(0); i++ ) {
	    double dif = fabs(epot(i)-epot2(i));
	    if( dif > eps ) {
		eps = dif;
		loc = i;
	    }
	}
	std::cout << "Maximum difference between comp. solutions: " << eps << " at " << loc << "\n";
    }
}


