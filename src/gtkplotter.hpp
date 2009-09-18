/*! \file gtkplotter.hpp
 *  \brief Header file for simulation plotter
 */


#ifndef GTK_PLOTTER_HPP
#define GTK_PLOTTER_HPP 1


#include <cairo.h>
#include <string>
#include <list>
#include <gtk/gtk.h>

#include "geometry.hpp"
#include "scalarfield.hpp"
#include "efield.hpp"
#include "vectorfield.hpp"
#include "particledatabase.hpp"
#include "particlediagplot.hpp"
#include "fielddiagplot.hpp"


/*! \page p_gtk_plotter GTK %Plotter
 *
 *  The interactive plotter is started by GTKPlotter. This object has
 *  to be created and given the data to plot with functions
 *  set_geometry(), set_epot(), etc. Then the a geometry plot can be
 *  started by calling new_geometry_plot_window() and run(). The
 *  simulation progress is halted while any plotter windows exist. See
 *  the following example:
 *
\code
    Geometry geom( MODE_2D, Int3D( 41, 41, 1 ), Vec3D( 0, 0, 0 ), 0.0025 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_DIRICHLET,  0.0) );
    geom.build_mesh();

    EpotProblem problem;
    problem.construct( geom );
    BiCGSTABSolver solver;
    problem.set_solver( solver );
    problem.solve( epot, scharge );

    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.new_geometry_plot_window();
    plotter.run();
\endcode
 *
 *  See class GTKPlotter for more information.
 */


/*! \brief GTK %Plotter class
 *
 *  The %GTKPlotter is the main simulation plotter class for
 *  interactive diagnostics.
 *
 */
class GTKPlotter {

    static bool             _gtk_initialized;

    std::list<class GTKWindow *>  _windows;
    
    const Geometry         *_geom;
    const ScalarField      *_epot;
    const Efield           *_efield;
    const ScalarField      *_scharge;
    const VectorField      *_bfield;
    const ParticleDataBase *_pdb;

public:

    GTKPlotter( int *argc, char ***argv );

    ~GTKPlotter();

    void run();

    const Geometry *get_geometry( void ) const;
    const ScalarField *get_epot( void ) const;
    const Efield *get_efield( void ) const;
    const ScalarField *get_scharge( void ) const;
    const VectorField *get_bfield( void ) const;
    const ParticleDataBase *get_particledatabase( void ) const;

    void set_geometry( const Geometry *geom );
    void set_epot( const ScalarField *epot );
    void set_efield( const Efield *efield );
    void set_scharge( const ScalarField *scharge );
    void set_bfield( const VectorField *bfield );
    void set_particledatabase( const ParticleDataBase *pdb );    

    GTKWindow *new_geometry_plot_window( void );
    GTKWindow *new_particle_plot_window( coordinate_axis_e axis, double level, 
					 particle_diag_plot_type_e type,
					 trajectory_diagnostic_e diagx, 
					 trajectory_diagnostic_e diagy );
    GTKWindow *new_field_plot_window( size_t N, const Vec3D &x1, const Vec3D &x2,
				      const field_diag_type_e diag[2], 
				      const field_loc_type_e loc[2] );

    //GTKWindow *new_window( void );
    void delete_window( class GTKWindow *window );
};


#endif













