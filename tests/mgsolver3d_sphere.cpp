/*! \file mgsolver3d_sphere.cpp 
 *  \brief Test multi grid solver with a 3d problem made of two concentric spheres.
 *
 *  \test Test multi grid solver with a 3d problem made of two concentric spheres.
 */


#include <fstream>
#include <iomanip>
#include "epot_bicgstabsolver.hpp"
#include "epot_gssolver.hpp"
#include "epot_mgsolver.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "gtkplotter.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "ibsimutest.hpp"


using namespace std;


double phi_a = 0.0;
double phi_b = 10.0;
double r_a = 0.02;
double r_b = 0.07;
double A = (phi_a - phi_b)/(1.0/r_a - 1.0/r_b);
double B = phi_a + A/r_a;


bool solid1( double x, double y, double z )
{
    return( x*x + y*y + z*z <= 0.02*0.02 );
}


bool solid2( double x, double y, double z )
{
    return( x*x + y*y + z*z >= 0.07*0.07 );
}


double phi( double r )
{
    if( r <= 0.02 )
	return( phi_a );
    else if( r >= 0.07 )
	return( phi_b );
    return( A/r - B );
}


void compare_to_analytic( const std::string &fieldname, const EpotField &epot )
{
    double eps = 0.0;
    size_t loci = 0, locj = 0, lock = 0;
    for( uint32_t i = 0; i < epot.size(0); i++ ) {
	for( uint32_t j = 0; j < epot.size(1); j++ ) {
	    for( uint32_t k = 0; k < epot.size(2); k++ ) {
		double x = epot.origo(0)+i*epot.h();
		double y = epot.origo(1)+j*epot.h();
		double z = epot.origo(2)+k*epot.h();
		double r = sqrt( x*x + y*y + z*z );
		if( r > 0.021 && r < 0.07 ) {
		    double dif = fabs( epot(i,j,k) - phi(r) );
		    if( dif > eps ) {
			eps = dif;
			loci = i;
			locj = j;
			lock = k;
		    }
		}
	    }
	}
    }
    std::cout << "Maximum difference between " << fieldname << " and analytic is " << eps << " at " 
	      << loci << ", " << locj << ", " << lock << "\n";
}


void test( int argc, char **argv )
{
    double h = 0.001;
    int32_t size = (int32_t)ceil(0.08/h) + 1;
    Geometry geom( MODE_3D, Int3D(size,size,size), Vec3D(0,0,0), h );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 10.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET, 10.0) );
    geom.set_boundary( 5, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 6, Bound(BOUND_DIRICHLET, 10.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, 10.0) );
    geom.build_mesh();
    
    EpotField epot( geom );
    MeshScalarField scharge( geom );

    EpotMGSolver solver( geom );
    solver.set_neumann_order( 2 );
    solver.set_levels( 5 );
    //solver.set_npre( 5 );
    //solver.set_npost( 2 );
    //solver.set_mgcycmax( 1 );
    solver.solve( epot, scharge );

    compare_to_analytic( "epot", epot );

    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.set_scharge( &scharge );
    plotter.set_epot( &epot );
    plotter.new_geometry_plot_window();
    plotter.run();
}

