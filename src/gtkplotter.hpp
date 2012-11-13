/*! \file gtkplotter.hpp
 *  \brief GTK based plotters.
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * taneli.kalvas@jyu.fi.
 * 
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
 */

#ifndef GTKPLOTTER_HPP
#define GTKPLOTTER_HPP 1


#include <vector>
#include <gtk/gtk.h>
#include "geometry.hpp"
#include "meshscalarfield.hpp"
#include "epot_field.hpp"
#include "epot_efield.hpp"
#include "vectorfield.hpp"
#include "particledatabase.hpp"
#include "particlediagplot.hpp"
#include "fielddiagplot.hpp"
#include "gtkwindow.hpp"


/*! \brief GTK %Plotter class
 *
 *  The %GTKPlotter is the main simulation plotter class for
 *  interactive diagnostics.
 *
 */
class GTKPlotter {

    static bool              _gtk_initialized;

    std::vector<GTKWindow *> _windows;
    
    const Geometry          *_geom;
    const EpotField         *_epot;
    const EpotEfield        *_efield;
    const MeshScalarField   *_scharge;
    const MeshScalarField   *_tdens;
    const VectorField       *_bfield;
    const ParticleDataBase  *_pdb;

public:

    GTKPlotter( int *argc, char ***argv );

    ~GTKPlotter();

    void run();

    const Geometry *get_geometry( void ) const;
    const EpotField *get_epot( void ) const;
    const EpotEfield *get_efield( void ) const;
    const MeshScalarField *get_scharge( void ) const;
    const MeshScalarField *get_trajdens( void ) const;
    const VectorField *get_bfield( void ) const;
    const ParticleDataBase *get_particledatabase( void ) const;

    void set_geometry( const Geometry *geom );
    void set_epot( const EpotField *epot );
    void set_efield( const EpotEfield *efield );
    void set_scharge( const MeshScalarField *scharge );
    void set_trajdens( const MeshScalarField *tdens );
    void set_bfield( const VectorField *bfield );
    void set_particledatabase( const ParticleDataBase *pdb );    

    GTKWindow *new_geometry_3d_plot_window( void );
    GTKWindow *new_geometry_plot_window( void );
    GTKWindow *new_particle_plot_window( coordinate_axis_e axis, double level, 
					 particle_diag_plot_type_e type,
					 trajectory_diagnostic_e diagx, 
					 trajectory_diagnostic_e diagy );
    GTKWindow *new_field_plot_window( size_t N, const Vec3D &x1, const Vec3D &x2,
				      const field_diag_type_e diag[2], 
				      const field_loc_type_e loc[2] );

    void delete_window( GTKWindow *window );
};


#endif

