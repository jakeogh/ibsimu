/*! \file emittancedef3d.cpp 
 *  \test Test beam definition by emittance in 3d
 */


#include <fstream>
#include <iomanip>
#include <iostream>
#include "geometry.hpp"
#include "epot_problem.hpp"
#include "epot_efield.hpp"
#include "bicgstab_solver.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"
#include "verbose.hpp"
#include "error.hpp"


using namespace std;


void test( void )
{
    Geometry geom( MODE_3D, Int3D(21,21,21), Vec3D(0,0,0), 1e-3 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 3, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 5, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 6, Bound(BOUND_DIRICHLET, 0.0) );
    geom.build_mesh();

    EpotProblem p;
    p.construct( geom );

    ScalarField epot( geom );
    ScalarField scharge( geom );

    BiCGSTABSolver solver;
    p.set_solver( solver );
    p.solve( epot, scharge );

    EpotEfield efield( geom, epot );
    VectorField bfield;

    ParticleDataBase3D pdb;
    bool pmirror[6] = { false, false, true, false, true, false };
    pdb.set_mirror( pmirror );

    // KV distribution

    pdb.add_3d_KV_beam_with_emittance( 10000, 10.0e-3, 6.0, 20.0, 
				       1.1, 0.193, 3.23e-5,
				       -1.1, 0.193, 6.6e-5,
				       3.0e3, 0.0, 0.0, 0.0 );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    // Plot (y,y') emittance and check values
    ParticleDiagPlotter pplotter1( &geom, &pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				  DIAG_Y, DIAG_YP );
    pplotter1.plot_png( "emittancedef3d_kv_yyp.png" );
    Emittance emit1 = pplotter1.calculate_emittance();
    if( fabs( 4.0*emit1.epsilon() - 3.23e-5 ) / 3.23e-5 > 0.01 )
	throw( Error( ERROR_LOCATION, "rms emittance does not match theory" ) );
    if( fabs( emit1.alpha() - 1.1 ) / 1.1 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted alpha does not match theory" ) );
    if( fabs( emit1.beta() - 0.193 ) / 0.193 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted beta does not match theory" ) );

    // Plot (z,z') emittance and check values
    ParticleDiagPlotter pplotter2( &geom, &pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				  DIAG_Z, DIAG_ZP );
    pplotter2.plot_png( "emittancedef3d_kv_zzp.png" );
    Emittance emit2 = pplotter2.calculate_emittance();
    if( fabs( 4.0*emit2.epsilon() - 6.6e-5 ) / 6.6e-5 > 0.01 )
	throw( Error( ERROR_LOCATION, "rms emittance does not match theory" ) );
    if( fabs( emit2.alpha() + 1.1 ) / 1.1 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted alpha does not match theory" ) );
    if( fabs( emit2.beta() - 0.193 ) / 0.193 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted beta does not match theory" ) );

    // Plot profile
    ParticleDiagPlotter pplotter3( &geom, &pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				  DIAG_Y, DIAG_Z );
    pplotter3.plot_png( "emittancedef3d_kv_yz.png" );


    // Gaussian distribution

    pdb.clear();
    pdb.add_3d_gaussian_beam_with_emittance( 100000, 10.0e-3, 6.0, 20.0, 
					     1.1, 0.193, 3.23e-5,
					     -1.1, 0.193, 6.6e-5,
					     3.0e3, 0.0, 0.0, 0.0 );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    // Plot (y,y') emittance and check values
    ParticleDiagPlotter pplotter4( &geom, &pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_Y, DIAG_YP );
    pplotter4.plot_png( "emittancedef3d_gaussian_yyp.png" );
    Emittance emit4 = pplotter4.calculate_emittance();
    if( fabs( emit4.epsilon() - 3.23e-5 ) / 3.23e-5 > 0.01 )
	throw( Error( ERROR_LOCATION, "rms emittance does not match theory" ) );
    if( fabs( emit4.alpha() - 1.1 ) / 1.1 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted alpha does not match theory" ) );
    if( fabs( emit4.beta() - 0.193 ) / 0.193 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted beta does not match theory" ) );

    // Plot (z,z') emittance and check value
    ParticleDiagPlotter pplotter5( &geom, &pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_Z, DIAG_ZP );
    pplotter5.plot_png( "emittancedef3d_gaussian_zzp.png" );
    Emittance emit5 = pplotter5.calculate_emittance();
    if( fabs( emit5.epsilon() - 6.6e-5 ) / 6.6e-5 > 0.01 )
	throw( Error( ERROR_LOCATION, "rms emittance does not match theory" ) );
    if( fabs( emit5.alpha() + 1.1 ) / 1.1 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted alpha does not match theory" ) );
    if( fabs( emit5.beta() - 0.193 ) / 0.193 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted beta does not match theory" ) );

    // Plot profile
    ParticleDiagPlotter pplotter6( &geom, &pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				  DIAG_Y, DIAG_Z );
    pplotter6.plot_png( "emittancedef3d_gaussian_yz.png" );
    
}


int main( void )
{
    try {
	verbose_output = 0;
	test();
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}

