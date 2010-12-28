/*! \file Test plotting to files png, eps, pdf, svg.
 *  \test Test plotting to files png, eps, pdf, svg.
 */


#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cairo.h>
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "vectorfield.hpp"
#include "particledatabase.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "geomplotter.hpp"


using namespace std;


bool solid1( double x, double y, double z )
{
    return( x <= 0.00187 && y >= 0.0015 && y >= 2.28*x - 0.0010 );
}


bool solid2( double x, double y, double z )
{
    return( x >= 0.0095 && y >= 0.0023333 && y >= 0.01283 - x );
}


void test( int *argc, char ***argv )
{    
    // 12x7 mm geometry with 0.05 mm mesh size
    Geometry geom( MODE_2D, Int3D(121,71,1), Vec3D(0,0,0), 0.0001 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0 ) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, -8.0e3) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0)  );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, -8.0e3) );
    geom.build_mesh();

    EpotProblem p;
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );

    VectorField bfield;
    EpotEfield efield( geom, epot );
    field_extrpl_e efldextrpl[6] = {FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				    FIELD_MIRROR,      FIELD_EXTRAPOLATE,
				    FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase2D pdb;
    pdb.set_thread_count( 4 );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );

    p.solve( epot, scharge );
    pdb.clear();
    pdb.add_2d_beam_with_energy( 5000, 600.0, 1.0, 1.0, 
				 500.0, 0.0, 0.5, 
				 0.0, 0.0, 
				 0.0, 0.0005 );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    GeomPlotter gplotter( &geom );
    gplotter.set_scharge( &scharge );
    //gplotter.set_scharge_field( true );
    gplotter.set_epot( &epot );
    gplotter.set_particle_database( &pdb );
    gplotter.set_particle_div( 100 );
    gplotter.set_font_size( 15 );

#ifdef CAIRO_HAS_PNG_FUNCTIONS
    gplotter.plot_png( "fileplot.png" );
#endif

#ifdef CAIRO_HAS_SVG_SURFACE
    gplotter.plot_svg( "fileplot.svg" );
#endif

#ifdef CAIRO_HAS_PS_SURFACE
    gplotter.plot_eps( "fileplot.eps" );
#endif

#ifdef CAIRO_HAS_PDF_SURFACE
    gplotter.plot_pdf( "fileplot.pdf" );
#endif

}


int main( int argc, char **argv )
{
    try {
	ibsimu.set_verbose_output( 0 );
	test( &argc, &argv );
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}



