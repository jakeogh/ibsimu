/*! \file mgsolver1d_dirichlet.cpp 
 *  \brief Test multi grid solver with a 1d problem with dirichlet boundaries
 *
 *  \test Test multi grid solver with a 1d problem with dirichlet boundaries
 */


#include <fstream>
#include <iomanip>
#include "epot_gssolver.hpp"
#include "epot_mgsolver.hpp"
#include "geometry.hpp"
#include "ibsimu.hpp"
#include "error.hpp"


using namespace std;


void test( int argc, char **argv )
{
    bool err = false;
    //Geometry geom( MODE_1D, Int3D(3,1,1), Vec3D(0,0,0), 0.02 ); // level 5
    //Geometry geom( MODE_1D, Int3D(5,1,1), Vec3D(0,0,0), 0.01 ); // level 4
    Geometry geom( MODE_1D, Int3D(9,1,1), Vec3D(0,0,0), 0.005 ); // level 3
    //Geometry geom( MODE_1D, Int3D(17,1,1), Vec3D(0,0,0), 0.0025 ); // level2
    //Geometry geom( MODE_1D, Int3D(33,1,1), Vec3D(0,0,0), 0.00125 ); // level 1
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 1000.0) );
    geom.build_mesh();

    MeshScalarField scharge( geom );
    for( size_t i = 0; i < scharge.size(0); i++ ) {
	scharge(i) = 1.0e-4;
    }

    EpotMGSolver solver( geom );
    MeshScalarField epot( geom );
    solver.set_levels( 2 );
    solver.set_gamma( 1 );
    solver.solve( epot, scharge );
    
    EpotGSSolver solver2( geom );
    MeshScalarField epot2( geom );
    solver2.solve( epot2, scharge );

    // Compare solutions
    std::cout << "epot2:\n";
    double eps = 0.0;
    size_t loc = 0;
    for( size_t i = 0; i < geom.size(0); i++ ) {
	std::cout << std::setw(8) << epot2(i) << " ";
	double dif = fabs(epot(i)-epot2(i));
	if( dif > eps ) {
	    eps = dif;
	    loc = i;
	}
    }
    std::cout << "\n";

    std::cout << "Maximum difference of " << eps << " at " << loc << "\n";

    if( err ) {
	std::cout << "Error: solved potentials differ from each other\n";
	exit( 1 );
    }

}


