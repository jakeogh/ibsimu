/*! \file gtkplotter.cpp
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

#include <sstream>
#include "ibsimu.hpp"
#include "gtkplotter.hpp"
#include "gtkwindow.hpp"
#include "gtkgeomwindow.hpp"
#include "gtkparticlediagwindow.hpp"
#include "gtkfielddiagwindow.hpp"


/* ************************************************ *
 * GTK PLOTTER
 * ************************************************ */



bool GTKPlotter::_gtk_initialized = false;


GTKPlotter::GTKPlotter( int *argc, char ***argv )
    : _geom(NULL), _epot(NULL), _scharge(NULL), _tdens(NULL), 
      _efield(NULL), _bfield(NULL), _pdb(NULL)
{
    if( !_gtk_initialized ) {
	if( gtk_init_check( argc, argv ) == FALSE )
	    throw( (int)1 );
	_gtk_initialized = true;
    }
}


GTKPlotter::~GTKPlotter()
{
    
}


void GTKPlotter::run()
{
    if( ibsimu.get_verbose_output() )
	std::cout << "Running GTKPlotter\n";

    gtk_main();

    if( ibsimu.get_verbose_output() )
	std::cout << "  Done\n";
}


GTKWindow *GTKPlotter::new_geometry_plot_window( void )
{
    GTKWindow *window = new GTKGeomWindow( this, _geom, _epot, _scharge, _tdens, _bfield, _pdb );
    _windows.push_back( window );

    return( window );
}

GTKWindow *GTKPlotter::new_particle_plot_window( coordinate_axis_e axis, double level, 
						 particle_diag_plot_type_e type,
						 trajectory_diagnostic_e diagx, 
						 trajectory_diagnostic_e diagy )
{
    GTKWindow *window = new GTKParticleDiagWindow( this, _pdb, _geom, axis, level, type, diagx, diagy );
    _windows.push_back( window );

    return( window );
}


GTKWindow *GTKPlotter::new_field_plot_window( size_t N, const Vec3D &x1, const Vec3D &x2,
					      const field_diag_type_e diag[2], 
					      const field_loc_type_e loc[2] )
{
    GTKWindow *window = new GTKFieldDiagWindow( this, _geom, N, x1, x2, diag, loc );
    _windows.push_back( window );

    return( window );
}



void GTKPlotter::delete_window( GTKWindow *window )
{
    std::list<GTKWindow *>::iterator it;

    //std::cout << "Delete window\n";

    for( it = _windows.begin(); it != _windows.end(); it++ ) {
	if( *it == window ) {
	    delete *it;
	    _windows.erase( it );
	    break;
	}
    }

    if( _windows.size() == 0 ) {
	//std::cout << "Last window deleted\n";
	gtk_main_quit();
    }
}


const Geometry *GTKPlotter::get_geometry( void ) const
{
    return( _geom );
}

const ScalarField *GTKPlotter::get_epot( void ) const
{
    return( _epot );
}


const VectorField *GTKPlotter::get_efield( void ) const
{
    return( _efield );
}


const ScalarField *GTKPlotter::get_scharge( void ) const
{
    return( _scharge );
}


const ScalarField *GTKPlotter::get_trajdens( void ) const
{
    return( _tdens );
}


const VectorField *GTKPlotter::get_bfield( void ) const
{
    return( _bfield );
}


const ParticleDataBase *GTKPlotter::get_particledatabase( void ) const
{
    return( _pdb );
}


void GTKPlotter::set_geometry( const Geometry *geom )
{
    _geom = geom;
}


void GTKPlotter::set_epot( const ScalarField *epot )
{
    _epot = epot;
}

void GTKPlotter::set_efield( const VectorField *efield )
{
    _efield = efield;
}

void GTKPlotter::set_scharge( const ScalarField *scharge )
{
    _scharge = scharge;
}


void GTKPlotter::set_trajdens( const ScalarField *tdens )
{
    _tdens = tdens;
}


void GTKPlotter::set_bfield( const VectorField *bfield )
{
    _bfield = bfield;
}


void GTKPlotter::set_particledatabase( const ParticleDataBase *pdb )
{
    _pdb = pdb;
}






