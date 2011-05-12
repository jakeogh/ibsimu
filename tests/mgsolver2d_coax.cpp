/*! \file mssolver2d_coax.cpp 
 *  \brief Test multi grid solver with a 2d problem made of two concentric cylinders.
 *
 *  \test Test multi grid solver with a 2d problem made of two concentric cylinders.
 */


#include <fstream>
#include <iomanip>
#include "epot_bicgstabsolver.hpp"
#include "epot_gssolver.hpp"
#include "epot_mgsolver.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_field.hpp"
#include "epot_efield.hpp"
#include "gtkplotter.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "ibsimutest.hpp"


using namespace std;


const double r1 = 0.021;
const double r2 = 0.07;
const double V1 = 0.0;
const double V2 = 10.0;
const double A = (V1-V2)/log(r1/r2);
const double B = V1-A*log(r1);


bool solid1( double x, double y, double z )
{
    return( x*x+y*y <= r1*r1 );
}


bool solid2( double x, double y, double z )
{
    return( x*x+y*y >= r2*r2 );
}


double phi( double r )
{
    return( A*log(r)+B );
}


void compare_to_analytic( const std::string &fieldname, const EpotField &epot )
{
    double eps = 0.0;
    size_t loci = 0, locj = 0;
    for( uint32_t i = 0; i < epot.size(0); i++ ) {
	for( uint32_t j = 0; j < epot.size(1); j++ ) {
	    double x = i*epot.h();
	    double y = j*epot.h();
	    double r = sqrt(x*x + y*y);
	    if( r > 0.021 && r < 0.07 ) {
	        double dif = fabs( epot(i,j) - phi(r) );
		if( dif > eps ) {
		    eps = dif;
		    loci = i;
		    locj = j;
		}
	    }
	}
    }
    std::cout << "Maximum difference between " << fieldname << " and analytic is " << eps << " at " << loci << ", " << locj << "\n";
}


void test( int argc, char **argv )
{
    //Geometry geom( MODE_2D, Int3D(6,6,1), Vec3D(0,0,0), 0.016 );        // one level possible
    //Geometry geom( MODE_2D, Int3D(11,11,1), Vec3D(0,0,0), 0.008 );      // two levels possible
    //Geometry geom( MODE_2D, Int3D(21,21,1), Vec3D(0,0,0), 0.004 );      // three levels possible
    //Geometry geom( MODE_2D, Int3D(41,41,1), Vec3D(0,0,0), 0.002 );      // four levels possible
    //Geometry geom( MODE_2D, Int3D(81,81,1), Vec3D(0,0,0), 0.001 );      // five levels possible
    //Geometry geom( MODE_2D, Int3D(161,161,1), Vec3D(0,0,0), 0.0005 );   // six levels possible
    //Geometry geom( MODE_2D, Int3D(321,321,1), Vec3D(0,0,0), 0.00025 );  // seven levels possible
    Geometry geom( MODE_2D, Int3D(641,641,1), Vec3D(0,0,0), 0.000125 );   // 8 levels possible
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,   V2) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,   V2) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,   V1) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET,   V2) );
    geom.build_mesh();
    //geom.debug_print( std::cout );

    MeshScalarField scharge( geom );

    EpotField epot( geom );
    EpotMGSolver solver( geom );
    solver.set_levels( 4 );
    solver.set_npre( 5 );
    solver.set_npost( 5 );
    solver.set_neumann_order( 2 );
    solver.solve( epot, scharge );

    ifstream is( "mgsolver2d_coax.dat" );
    EpotField epot2( is, geom );
    is.close();

    /*
    //EpotField epot2( geom );
    EpotGSSolver solver2( geom );
    solver2.set_w( 1.00 );
    solver2.set_eps( 1.0e-8 );
    solver2.set_neumann_order( 2 );
    solver2.solve( epot2, scharge );
    */

    // Save GS solution
    //ofstream os( "mgsolver2d_coax.dat" );
    //epot2.save( os );
    //os.close();

    /*
    EpotField epot3( geom );
    EpotBiCGSTABSolver solver3( geom );
    solver3.set_neumann_order( 2 );
    solver3.solve( epot3, scharge );
    */

    // Compare numeric solutions
    double eps = 0.0;
    size_t loci = 0, locj = 0;
    for( size_t j = 0; j < geom.size(1); j++ ) {
	for( size_t i = 0; i < geom.size(0); i++ ) {
	    double dif = fabs(epot(i,j)-epot2(i,j));
	    if( dif > eps ) {
		eps = dif;
		loci = i;
		locj = j;
	    }
	}
    }
    std::cout << "Maximum difference between epot and epot2 is " << eps << " at " << loci << ", " << locj << "\n";

    compare_to_analytic( "epot", epot );
    compare_to_analytic( "epot2", epot2 );

    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.new_geometry_plot_window();
    plotter.run();
}


