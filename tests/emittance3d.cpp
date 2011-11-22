/*! \file emittance3d.cpp 
 *  \brief Test emittance of 3d beam against analytical value
 * 
 *  \test Test emittance of 3d beam against analytical value
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
    Geometry geom( MODE_3D, Int3D(21,21,21), Vec3D(0,-1e-2,-1e-2), 1e-3 );
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
    MeshVectorField bfield;

    // Make a circular beam of Neon (m=20) 6+ with J = 45 A/m2 and radius 5 mm.
    // Starting energy is 3 keV and transverse temperature 1.0 eV.
    // The emittance should is known analytically and is 0.5*sqrt(kT/m)*r/v_z = 
    // 3.22748612e-5 pi*m*rad 
    ParticleDataBase3D pdb;
    pdb.add_cylindrical_beam_with_energy( 20000, 45.0, 6.0, 20.0, 
					  3.0e3, 0.0, 1.0, 
					  Vec3D(0,0,0),
					  Vec3D(0,1,0),
					  Vec3D(0,0,1),
					  0.005 );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    // Plot (y,y') emittance and check value
    ParticleDiagPlotter pplotter( geom, pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				  DIAG_Y, DIAG_YP );
    pplotter.plot_png( "emittance3d_emit.png" );
    Emittance emit = pplotter.calculate_emittance();
    if( fabs( emit.epsilon() - 3.22748612e-5 ) / 3.22748612e-5 > 0.01 )
	throw( Error( ERROR_LOCATION, "rms emittance does not match theory" ) );

    // Plot (y,z) profile
    ParticleDiagPlotter pplotter2( geom, pdb, AXIS_X, 1.0e-6, PARTICLE_DIAG_PLOT_HISTO2D,
				   DIAG_Y, DIAG_Z );
    pplotter2.plot_png( "emittance3d_profile.png" );
}



