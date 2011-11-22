/*! \file emittancedef2d.cpp 
 *  \brief Test beam definition by emittance in 2d
 *
 *  \test Test beam definition by emittance in 2d
 */


#include <fstream>
#include <iomanip>
#include <iostream>
#include "geometry.hpp"
#include "meshvectorfield.hpp"
#include "epot_problem.hpp"
#include "epot_efield.hpp"
#include "bicgstab_solver.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"
#include "ibsimu.hpp"
#include "error.hpp"


using namespace std;


void test( int argc, char **argv )
{
    Geometry geom( MODE_2D, Int3D(21,21,1), Vec3D(0,0,0), 1e-3 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 3, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET, 0.0) );
    geom.build_mesh();

    EpotProblem p;
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    EpotEfield efield( geom, epot );
    MeshVectorField bfield;

    ParticleDataBase2D pdb;
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );

    // KV distribution

    pdb.add_2d_KV_beam_with_emittance( 10000, 10.0e-3, 6.0, 20.0, 
				       1.1, 0.193, 3.23e-5,
				       3.0e3, 0.0, 0.0 );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    // Plot (y,y') emittance and check value
    ParticleDiagPlotter pplotter( geom, pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				  DIAG_Y, DIAG_YP );
    pplotter.set_ranges( -0.006, -0.0601, 0.006, 0.0601 );
    pplotter.plot_png( "emittancedef2d_kv.png" );
    Emittance emit = pplotter.calculate_emittance();
    if( fabs( emit.epsilon() - 3.23e-5 ) / 3.23e-5 > 0.01 )
	throw( Error( ERROR_LOCATION, "rms emittance does not match theory" ) );
    if( fabs( emit.alpha() - 1.1 ) / 1.1 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted alpha does not match theory" ) );
    if( fabs( emit.beta() - 0.193 ) / 0.193 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted beta does not match theory" ) );


    // Gaussian distribution

    pdb.clear();
    pdb.add_2d_gaussian_beam_with_emittance( 10000, 10.0e-3, 6.0, 20.0, 
					     1.1, 0.193, 3.23e-5,
					     3.0e3, 0.0, 0.0 );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    // Plot (y,y') emittance and check value
    ParticleDiagPlotter pplotter2( geom, pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_Y, DIAG_YP );
    //pplotter2.set_ranges( -0.006, -0.0601, 0.006, 0.0601 );
    pplotter2.plot_png( "emittancedef2d_gaussian.png" );
    Emittance emit2 = pplotter.calculate_emittance();
    if( fabs( emit2.epsilon() - 3.23e-5 ) / 3.23e-5 > 0.01 )
	throw( Error( ERROR_LOCATION, "rms emittance does not match theory" ) );
    if( fabs( emit2.alpha() - 1.1 ) / 1.1 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted alpha does not match theory" ) );
    if( fabs( emit2.beta() - 0.193 ) / 0.193 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted beta does not match theory" ) );

    
}



