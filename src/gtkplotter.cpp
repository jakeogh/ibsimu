#include <sstream>
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
    : _geom(NULL), _epot(NULL), _efield(NULL), 
      _scharge(NULL), _bfield(NULL), _pdb(NULL)
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
    gtk_main();
}


GTKWindow *GTKPlotter::new_geometry_plot_window( void )
{
    GTKWindow *window = new GTKGeomWindow( this, _geom, _epot, _scharge, _bfield, _pdb );
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


const Efield *GTKPlotter::get_efield( void ) const
{
    return( _efield );
}


const ScalarField *GTKPlotter::get_scharge( void ) const
{
    return( _scharge );
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

void GTKPlotter::set_efield( const Efield *efield )
{
    _efield = efield;
}

void GTKPlotter::set_scharge( const ScalarField *scharge )
{
    _scharge = scharge;
}


void GTKPlotter::set_bfield( const VectorField *bfield )
{
    _bfield = bfield;
}


void GTKPlotter::set_particledatabase( const ParticleDataBase *pdb )
{
    _pdb = pdb;
}













