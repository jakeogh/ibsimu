/*! \file energy.cpp 
 *  \brief Test energy distribution with a beam in 2d/cyl system.
 *
 *  \test Test energy distribution with a beam in 2d/cyl system.
 */


#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_bicgstabsolver.hpp"
#include "epot_field.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"
#include "ibsimu.hpp"
#include "ibsimutest.hpp"
#include "error.hpp"
#include "gtkplotter.hpp"


using namespace std;


void parallel( int argc, char **argv )
{    
    // 10x10 cm geometry with 0.25 cm mesh
    Geometry geom( MODE_2D, Int3D(41,41,1), Vec3D(0,0,0), 0.0025 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,  0.0) );
    geom.build_mesh();

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    EpotEfield efield( epot );
    MeshVectorField bfield;

    ParticleDataBase2D pdb;
    pdb.set_thread_count( 1 );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.add_2d_beam_with_energy( 10000, 50.0, 1.0, 1.0, 
				 3.0e3, 1.0, 0.0, 
				 0.0, 0.0, 
				 0.0, 0.009 );
    pdb.set_accuracy( 1.0e-9, 1.0e-9 );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    ParticleDiagPlotter diag( geom, pdb, AXIS_X, 0.0, PARTICLE_DIAG_PLOT_HISTO1D, 
			      DIAG_VX );
    diag.set_histogram_n( 100 );
    diag.plot_png( "energy_pvel.png" );
    diag.export_data( "energy_pvel.dat" );

    ParticleDiagPlotter diag2( geom, pdb, AXIS_X, 0.0, PARTICLE_DIAG_PLOT_HISTO1D, 
			      DIAG_EK );
    diag.set_histogram_n( 100 );
    diag2.plot_png( "energy_pen.png" );
    diag2.export_data( "energy_pen.dat" );

    if( false ) {
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.set_scharge( &scharge );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
    }
}


void transverse( int argc, char **argv )
{    
    // 10x10 cm geometry with 0.25 cm mesh
    Geometry geom( MODE_2D, Int3D(41,41,1), Vec3D(0,0,0), 0.0025 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,  0.0) );
    geom.build_mesh();

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    EpotEfield efield( epot );
    MeshVectorField bfield;

    ParticleDataBase2D pdb;
    pdb.set_thread_count( 1 );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.add_2d_beam_with_energy( 10000, 50.0, 1.0, 1.0, 
				 3.0e3, 0.0, 1.0, 
				 0.0, 0.0, 
				 0.0, 0.009 );
    pdb.set_accuracy( 1.0e-9, 1.0e-9 );
    pdb.iterate_trajectories( scharge, efield, bfield, geom );

    ParticleDiagPlotter diag( geom, pdb, AXIS_X, 0.0, PARTICLE_DIAG_PLOT_HISTO1D, 
			      DIAG_VY );
    diag.set_histogram_n( 100 );
    diag.plot_png( "energy_tvel.png" );
    diag.export_data( "energy_tvel.dat" );

    if( false ) {
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.set_scharge( &scharge );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
    }
}


void total( int argc, char **argv )
{    
    double h = 1.0e-4;
    double sizereq[3] = { 200.0e-3,
                           15.0e-3, 
                            0.0e-3 };
    Int3D meshsize( (int)floor(sizereq[0]/h)+1,
                    (int)floor(sizereq[1]/h)+1,
                    (int)floor(sizereq[2]/h)+1 );
    Vec3D origo( 0.0, 0.0, 0.0 );
    Geometry geom( MODE_CYL, meshsize, origo, h );

    geom.set_boundary( 1, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,  0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,  0.0) );
    geom.build_mesh();

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    EpotEfield efield( epot );
    MeshVectorField bfield;

    EpotBiCGSTABSolver solver( geom );

    ParticleDataBaseCyl pdb;
    pdb.set_thread_count( 1 );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_accuracy( 1.0e-9, 1.0e-9 );

    const double m = 1.0;
    const double q = 1.0;
    const double r0 = 5.0e-3;
    const double I = 5.0e-3;
    const double J = I/(M_PI*r0*r0);
    const double Etot = 20.0e3;
    for( uint32_t a = 0; a < 10; a++ ) {
	
	solver.solve( epot, scharge );
        efield.recalculate();

	pdb.clear();
	pdb.add_2d_beam_with_total_energy( 10000, J, q, m, 
					   Etot, epot, 0.0, 0.0, 
					   0.0, 0.0, 
					   0.0, r0 );
	pdb.iterate_trajectories( scharge, efield, bfield, geom );
    }

    // Calculate total energy
    TrajectoryDiagnosticData tdata;
    std::vector<trajectory_diagnostic_e> diagnostics;
    diagnostics.push_back( DIAG_X );
    diagnostics.push_back( DIAG_R );
    diagnostics.push_back( DIAG_VX );
    diagnostics.push_back( DIAG_VR );
    diagnostics.push_back( DIAG_VTHETA );
    pdb.trajectories_at_plane( tdata, AXIS_X, 0.0, diagnostics );
    //std::vector<double> datar;
    //std::vector<double> dataEtot;
    ofstream dout( "energy_etot.dat" );
    dout << "# "
	 << setw(10) << "r (m)" << " "
	 << setw(12) << "Ekin (eV)" << " "
	 << setw(12) << "Epot (eV)" << "\n";
    double errmax = 0.0;
    for( uint32_t a = 0; a < tdata.traj_size(); a++ ) {
	double r = tdata(a,1);
	double v2 = tdata(a,2)*tdata(a,2) + 
	    tdata(a,3)*tdata(a,3) +
	    tdata(a,4)*tdata(a,4);
	double Ekin = 0.5*MASS_U*m*v2;
	Vec3D loc( tdata(a,0), tdata(a,1), 0.0 );
	double U = epot( loc );
	double Epot = CHARGE_E*q*U;
	double err = fabs(20e3 - (Ekin+Epot)/CHARGE_E);
	if( err > errmax )
	    errmax = err;
	    
	dout << setw(12) << r << " "
	     << setw(12) << Ekin/CHARGE_E << " "
	     << setw(12) << Epot/CHARGE_E << "\n";
    }

    if( errmax > 1.0 )
	throw( ErrorTest( ERROR_LOCATION, "Kinetic error differs from expected by " + to_string(errmax) ) );
    
    if( false ) {
	GTKPlotter plotter( &argc, &argv );
	plotter.set_geometry( &geom );
	plotter.set_epot( &epot );
	plotter.set_scharge( &scharge );
	plotter.set_particledatabase( &pdb );
	plotter.new_geometry_plot_window();
	plotter.run();
    }
}


void test( int argc, char **argv )
{
    parallel( argc, argv );
    transverse( argc, argv );
    total( argc, argv );
}
