/*! \file emittancedefcyl.cpp 
 *  \test Test beam definition by emittance in cyl
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
#include "ibsimu.hpp"
#include "error.hpp"


using namespace std;


void test( void )
{
    Geometry geom( MODE_CYL, Int3D(21,21,1), Vec3D(0,0,0), 1e-3 );
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
    VectorField bfield;

    ParticleDataBaseCyl pdb;
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );

    pdb.clear();
    pdb.add_2d_gaussian_beam_with_emittance( 10000, 10.0e-3, 6.0, 20.0, 
					     1.1, 0.193, 3.23e-5,
					     3.0e3, 0.0 );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    // Plot (r,r') emittance
    ParticleDiagPlotter pplotter1( &geom, &pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_R, DIAG_RP );
    pplotter1.plot_png( "emittancedefcyl_gaussian_rrp.png" );

    // Plot (r,a') emittance
    ParticleDiagPlotter pplotter2( &geom, &pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_R, DIAG_AP );
    pplotter2.plot_png( "emittancedefcyl_gaussian_rap.png" );

    // Plot (r',a') emittance
    ParticleDiagPlotter pplotter3( &geom, &pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_RP, DIAG_AP );
    pplotter3.plot_png( "emittancedefcyl_gaussian_rpap.png" );

    // Plot (y,y') emittance and check value
    ParticleDiagPlotter pplotter4( &geom, &pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_Y, DIAG_YP );
    pplotter4.plot_png( "emittancedefcyl_gaussian_yyp.png" );
    Emittance emit4 = pplotter4.calculate_emittance();
    if( fabs( emit4.epsilon() - 3.23e-5 ) / 3.23e-5 > 0.01 )
	throw( Error( ERROR_LOCATION, "rms emittance does not match theory" ) );
    if( fabs( emit4.alpha() - 1.1 ) / 1.1 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted alpha does not match theory" ) );
    if( fabs( emit4.beta() - 0.193 ) / 0.193 > 0.01 )
	throw( Error( ERROR_LOCATION, "fitted beta does not match theory" ) );

    
}


int main( void )
{
    try {
	ibsimu.set_verbose_output( 0 );
	test();
    } catch ( Error e ) {
	cout << "Error in " << e._loc._file << ":" << e._loc._line 
	     << " in " << e._loc._func << "(): " << e._error_str << "\n";
	exit( 1 );
    }

    return( 0 );
}



