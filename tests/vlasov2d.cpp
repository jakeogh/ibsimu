/*! \file Test with vlasov iteration in 2d electrode configuration.
 *  \brief Test with vlasov iteration in 2d electrode configuration.
 *
 *  \test Test with vlasov iteration in 2d electrode configuration.
 */


#include <iostream>
#include <iomanip>
#include "epot_gssolver.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_field.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "particledatabase.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "ibsimu.hpp"
#include "error.hpp"


using namespace std;


bool solid1( double x, double y, double z )
{
    return( x <= 0.02 && y >= 0.018 );
}


bool solid2( double x, double y, double z )
{
    return( x >= 0.03 && x <= 0.04 && y >= 0.02 );
}


bool solid3( double x, double y, double z )
{
    return( x >= 0.06 && y >= 0.03 && y >= 0.07 - 0.5*x );
}


void test( int argc, char **argv )
{
    // 12x5 cm geometry with 0.05 cm mesh size
    Geometry geom( MODE_2D, Int3D(241,101,1), Vec3D(0,0,0), 0.0005 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    Solid *s3 = new FuncSolid( solid3 );
    geom.set_solid( 9, s3 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,  -3.0e3) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,  -1.0e3) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,     0.0  ) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,     0.0  ) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  -3.0e3) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, -14.0e3) );
    geom.set_boundary( 9, Bound(BOUND_DIRICHLET,  -1.0e3) );
    geom.build_mesh();

    EpotGSSolver solver( geom );
    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;
    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_MIRROR,      FIELD_EXTRAPOLATE,
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase2D pdb;
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );

    for( size_t i = 0; i < 10; i++ ) {
	solver.solve( epot, scharge );
	efield.recalculate();
	pdb.clear();
	pdb.add_2d_beam_with_energy( 1000, 50.0, 1.0, 1.0, 
				     3.0e3, 0.0, 0.0, 
				     0.0, 0.0, 
				     0.0, 0.012 );
	pdb.iterate_trajectories( scharge, efield, bfield, geom );

	/*
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.set_efield( &efield );
	plotter.set_scharge( &scharge );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
	*/
    }

    GeomPlotter geomplotter( geom );
    geomplotter.set_size( 1024, 768 );
    geomplotter.set_epot( &epot );
    std::vector<double> eqlines;
    eqlines.push_back( -5.0 );
    eqlines.push_back( 0.0 );
    eqlines.push_back( +5.0 );
    geomplotter.set_eqlines_manual( eqlines );
    geomplotter.set_particle_database( &pdb );
    geomplotter.plot_png( "vlasov2d.png" );
}


