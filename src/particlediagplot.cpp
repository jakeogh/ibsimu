/*! \file particlediagplot.cpp
 *  \brief Source code for particlediagplot.cpp
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

#include <limits>
#include <fstream>
#include "particlediagplot.hpp"
#include "ibsimu.hpp"


ParticleDiagPlot::ParticleDiagPlot( Frame *frame, const Geometry *geom, const ParticleDataBase *pdb, 
				    coordinate_axis_e axis, double level, 
				    particle_diag_plot_type_e type,
				    trajectory_diagnostic_e diagx, trajectory_diagnostic_e diagy )
    : _frame(frame), _geom(geom), _pdb(pdb), _axis(axis), _level(level), 
      _type(type), _diagx(diagx), _diagy(diagy), _diagz(DIAG_NONE),
      _pdb_it_no(-1), _update(true), _tdata(NULL), _histo(NULL), _emit(NULL),
      _scatter(NULL), _ellipse(NULL), _ellipse_enable(true), _colormap(NULL), _profile(NULL), 
      _histogram_n(50), _histogram_m(50), _interpolation(INTERPOLATION_CLOSEST), _dot_size(1.0)
{

}


ParticleDiagPlot::~ParticleDiagPlot()
{
    if( _tdata )
	delete _tdata;
    if( _histo )
	delete _histo;
    if( _emit )
	delete _emit;

    if( _scatter )
	delete _scatter;
    if( _colormap )
	delete _colormap;
    if( _ellipse )
	delete _ellipse;
    if( _profile )
	delete _profile;
}


/* Build data for particle diagnostics only if update needed (particle 
 * iteration redone or plot parameters changed)
 *
 * 1. Only in case of scatter plot return tdata
 * 2. Only in case of histogram plot return histogram
 * 3. Always when applicable, return emittance fit
 */
void ParticleDiagPlot::build_data( void )
{
    // Is update needed?
    if( _pdb->get_iteration_number() == _pdb_it_no && !_update )
	return;

    _pdb_it_no = _pdb->get_iteration_number();
    _update = false;

    // Free old data
    if( _tdata )
	delete _tdata;
    if( _histo )
	delete _histo;
    if( _emit )
	delete _emit;
    _tdata = NULL;
    _histo = NULL;
    _emit  = NULL;

    // Prepare diagnostic request
    bool emittanceconv = false;
    std::vector<trajectory_diagnostic_e> diagnostics;
    if( _type == PARTICLE_DIAG_PLOT_HISTO2D && _geom->geom_mode() == MODE_CYL && 
	( (_diagx == DIAG_Y && _diagy == DIAG_YP) || 
	  (_diagx == DIAG_Z && _diagy == DIAG_ZP)) ) {

	// Diagnostics for emittance conversion
	emittanceconv = true;
	diagnostics.push_back( DIAG_R );
	diagnostics.push_back( DIAG_RP );
	diagnostics.push_back( DIAG_AP );
	diagnostics.push_back( DIAG_CURR );

    } else {

	if( _diagx == DIAG_NONE )
	    throw( Error( ERROR_LOCATION, "no first diagnostic" ) );	
	diagnostics.push_back( _diagx );
	if( _diagy != DIAG_NONE )
	    diagnostics.push_back( _diagy );
	if( _diagz != DIAG_NONE )
	    diagnostics.push_back( _diagz );
	diagnostics.push_back( DIAG_CURR );
    }

    // Get diagnostic data
    _tdata = new TrajectoryDiagnosticData;
    _pdb->trajectories_at_plane( *_tdata, _axis, _level, diagnostics );

    // Do data mirroring. Limited to only one mirroring per
    // axis-direction, lower end dominates if both edges have
    // mirroring enabled. No mirroring done if making emittance
    // conversion plot.
    if( !emittanceconv ) {
	bool mirror[6];
	_pdb->get_mirror( mirror );
	// Mirror in x-direction
	if( mirror[0] ) {
	    _tdata->mirror( AXIS_X, _geom->origo(0) );
	} else if( mirror[1] ) {
	    _tdata->mirror( AXIS_X, _geom->max(0) );
	}
	// Mirror in y-direction, forced if cylindrical geometry
	if( _geom->geom_mode() == MODE_CYL || mirror[2] ) {
	    _tdata->mirror( AXIS_Y, _geom->origo(1) );
	} else if( mirror[3] ) {
	    _tdata->mirror( AXIS_Y, _geom->max(1) );
	}
	// Mirror in z-direction
	if( mirror[4] ) {
	    _tdata->mirror( AXIS_Z, _geom->origo(2) );
	} else if( mirror[5] ) {
	    _tdata->mirror( AXIS_Z, _geom->max(2) );
	}
    }

    if( emittanceconv ) {
	
	// Make emittance conversion histogram
	EmittanceConv *emit = new EmittanceConv( _histogram_n, _histogram_m,
						 (*_tdata)(0).data(), (*_tdata)(1).data(), 
						 (*_tdata)(2).data(), (*_tdata)(3).data() );

	// Copy and save emittance fit and histogram
	_histo = new Histogram2D( emit->histogram() );
	_emit  = emit;
	emit->free_histogram();

	// Delete invalid tdata
	delete _tdata;
	_tdata = NULL;
	
    } else if( _type == PARTICLE_DIAG_PLOT_SCATTER ) {

	// Build emittance fit if applicable
	if( (_diagx == DIAG_X && _diagy == DIAG_XP)  ||
	    (_diagx == DIAG_Y && _diagy == DIAG_YP)  ||
	    (_diagx == DIAG_R && _diagy == DIAG_RP)  ||
	    (_diagx == DIAG_Z && _diagy == DIAG_ZP) ) {
	    _emit = new Emittance( (*_tdata)(0).data(), 
				   (*_tdata)(1).data(), 
				   (*_tdata)(2).data() );
	}

    } else if( _type == PARTICLE_DIAG_PLOT_HISTO2D ) {

	// Make colormap
	Histogram2D *histo2d = new Histogram2D( _histogram_n, _histogram_m, (*_tdata)(0).data(), 
						(*_tdata)(1).data(), (*_tdata)(2).data() );
	_histo = histo2d;

	if( _geom->geom_mode() == MODE_CYL && _diagx == DIAG_R ) {
	    // Scale plot to have constant area per square for cylindrical geometry.
	    double dr = histo2d->nstep();
	    for( size_t i = 0; i < histo2d->n(); i++ ) {
		double r = fabs( histo2d->icoord( i ) );
		double w = M_PI*((r+0.5*dr)*(r+0.5*dr) - (r-0.5*dr)*(r-0.5*dr));
		if( r == 0.0 )
		    w = M_PI*dr*dr;
		for( size_t j = 0; j < histo2d->m(); j++ )
		    (*histo2d)(i,j) /= w;
	    }
	} else if( (_diagx == DIAG_X || _diagx == DIAG_Y || _diagx == DIAG_R || _diagx == DIAG_Z) && 
		   (_diagy == DIAG_X || _diagy == DIAG_Y || _diagy == DIAG_R || _diagy == DIAG_Z) ) {
	    // Profile plot: scale for A/m2 unit
	    double w = 1.0/(histo2d->nstep()*histo2d->mstep());
	    for( size_t i = 0; i < histo2d->n(); i++ )
		for( size_t j = 0; j < histo2d->m(); j++ )
		    (*histo2d)(i,j) *= w;
	}

	// Build emittance fit if applicable
	if( (_diagx == DIAG_X && _diagy == DIAG_XP)  ||
	    (_diagx == DIAG_Y && _diagy == DIAG_YP)  ||
	    (_diagx == DIAG_R && _diagy == DIAG_RP)  ||
	    (_diagx == DIAG_Z && _diagy == DIAG_ZP) ) {
	    _emit = new Emittance( (*_tdata)(0).data(), 
				   (*_tdata)(1).data(), 
				   (*_tdata)(2).data() );
	}

	// Delete unnecessary tdata
	delete _tdata;
	_tdata = NULL;

    } else if( _type == PARTICLE_DIAG_PLOT_HISTO1D ) {

	// Scale trajectory currents by 1/(2*pi*r)
	if( _geom->geom_mode() == MODE_CYL && _diagx == DIAG_R ) {
	    for( size_t a = 0; a < _tdata->traj_size(); a++ ) {
		if( (*_tdata)(a,0) != 0.0 )
		    (*_tdata)(a,1) /= (2.0*M_PI*fabs((*_tdata)(a,0)));
		else
		    (*_tdata)(a,1) = 0.0;
	    }
	}

	// Make XYGraph profile plot
	Histogram1D *histo1d = new Histogram1D( _histogram_n, (*_tdata)(0).data(), (*_tdata)(1).data() );
	_histo = histo1d;

	// Scale profile plot to have constant area per bin.
	if( _geom->geom_mode() == MODE_CYL && _diagx == DIAG_R ) {
	    double dr = histo1d->step();
	    for( size_t i = 0; i < histo1d->n(); i++ ) {
		(*histo1d)(i) /= dr;
	    }
	}

	// Delete unnecessary tdata
	delete _tdata;
	_tdata = NULL;

    }
}


const Emittance &ParticleDiagPlot::calculate_emittance( void )
{
    build_data();

    if( _emit )
	return( *_emit );

    throw( Error( ERROR_LOCATION, "diagnostic not emittance" ) );
}


void ParticleDiagPlot::export_data( const std::string &filename )
{
    build_data();

    std::ofstream fstr( filename.c_str() );
    if( ibsimu.get_verbose_output() )
	std::cout << "Exporting particle diagnostic data to \'" << filename << "\'\n";

    // Write header
    if( _type == PARTICLE_DIAG_PLOT_HISTO1D ) {
	fstr << "# ";
	fstr << std::setw(11) << trajectory_diagnostic_string_with_unit[_diagx];
	if( _geom->geom_mode() == MODE_CYL && _diagx == DIAG_R )
	    fstr << std::setw(14) << "Int (A/m^2)";
	else
	    fstr << std::setw(14) << "Int (a.u.)";
	fstr << "\n";
    } else {
	fstr << "# ";
	fstr << std::setw(11) << trajectory_diagnostic_string_with_unit[_diagx];
	fstr << std::setw(14) << trajectory_diagnostic_string_with_unit[_diagy];
	if( _geom->geom_mode() == MODE_CYL && _diagx == DIAG_R && _diagy == DIAG_RP)
	    fstr << std::setw(14) << "Int (A/m/rad)";
	else
	    fstr << std::setw(14) << "Int (a.u.)";
	fstr << "\n";
    }

    // Write data
    if( _type == PARTICLE_DIAG_PLOT_SCATTER ) {
	size_t N = _tdata->traj_size();
	for( size_t a = 0; a < N; a++ ) {
	    fstr << std::setw(13) << (*_tdata)(a,0) << " ";
	    fstr << std::setw(13) << (*_tdata)(a,1) << " ";
	    fstr << std::setw(13) << (*_tdata)(a,2) << "\n";
	}
    } else if( _type == PARTICLE_DIAG_PLOT_HISTO1D ) {
	Histogram1D *histo1d = dynamic_cast<Histogram1D *>( _histo );
	if( !histo1d )
	    throw( Error( ERROR_LOCATION, "dynamic cast error" ) );
	size_t N = histo1d->n();
	for( size_t a = 0; a < N; a++ ) {
	    fstr << std::setw(13) << histo1d->coord(a) << " ";
	    fstr << std::setw(13) << (*histo1d)(a) << "\n";
	}
    } else if( _type == PARTICLE_DIAG_PLOT_HISTO2D ) {
	Histogram2D *histo2d = dynamic_cast<Histogram2D *>( _histo );
	if( !histo2d )
	    throw( Error( ERROR_LOCATION, "dynamic cast error" ) );
	size_t N = histo2d->n();
	size_t M = histo2d->m();
	for( size_t b = 0; b < M; b++ ) {
	    for( size_t a = 0; a < N; a++ ) {
		fstr << std::setw(13) << histo2d->icoord(a) << " ";
		fstr << std::setw(13) << histo2d->jcoord(b) << " ";
		fstr << std::setw(13) << (*histo2d)(a,b) << "\n";
	    }
	    fstr << "\n";
	}
    }

    fstr.close();
}


void ParticleDiagPlot::build_plot( void )
{
    build_data();

    // Set axis labels
    if( _type == PARTICLE_DIAG_PLOT_HISTO1D ) {
	_frame->set_axis_label( PLOT_AXIS_X1, trajectory_diagnostic_string_with_unit[_diagx] );
	_frame->set_axis_label( PLOT_AXIS_Y1, "Intensity (a.u.)" );
    } else {
	_frame->set_axis_label( PLOT_AXIS_X1, trajectory_diagnostic_string_with_unit[_diagx] );
	_frame->set_axis_label( PLOT_AXIS_Y1, trajectory_diagnostic_string_with_unit[_diagy] );
    }

    // Clear old graphs
    _frame->clear_graphs();
    if( _scatter ) {
	delete _scatter;
	_scatter = NULL;
    }
    if( _colormap ) {
	delete _colormap; 
	_colormap = NULL;
    }
    if( _ellipse ) {
	delete _ellipse;
	_ellipse = NULL;
    }
    if( _profile ) {
	delete _profile;
	_profile = NULL;
    }

    if( _type == PARTICLE_DIAG_PLOT_SCATTER ) {

	// Make scatter
	_scatter = new XYGraph( (*_tdata)(0).data(), (*_tdata)(1).data() );
	_scatter->set_line_style( XYGRAPH_LINE_DISABLE );
	_scatter->set_point_style( XYGRAPH_POINT_CIRCLE, true, _dot_size );
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _scatter );

    } else if( _type == PARTICLE_DIAG_PLOT_HISTO2D ) {

	// Make colormap
	Histogram2D *histo2d = dynamic_cast<Histogram2D *>( _histo );
	if( !histo2d )
	    throw( Error( ERROR_LOCATION, "dynamic cast error" ) );
	double range[4];
	histo2d->get_range( range );
	double zmin, zmax;
	histo2d->get_bin_range( zmin, zmax );

	Palette palette;
	palette.clear();
	if( zmin == 0.0 && zmax >= 0.0 ) {
	    // Palette for positive beam
	    palette.push_back( Color(1,1,1), 0 );
	    palette.push_back( Color(1,1,0), 1 );
	    palette.push_back( Color(1,0,0), 2 );
	    palette.push_back( Color(0,0,0), 3 );
	} else if( zmax == 0.0 && zmin <= 0.0 ) {
	    // Palette for negative beam
	    palette.push_back( Color(1,1,1), 3 );
	    palette.push_back( Color(1,1,0), 2 );
	    palette.push_back( Color(1,0,0), 1 );
	    palette.push_back( Color(0,0,0), 0 );
	} else {
	    // Palette for positive and negative beam
	    palette.push_back( Color(0,0,0), zmin );
	    palette.push_back( Color(0,0,1), 0.67*zmin );
	    palette.push_back( Color(0,1,1), 0.33*zmin );
	    palette.push_back( Color(1,1,1), 0 );
	    palette.push_back( Color(1,1,0), 0.33*zmax );
	    palette.push_back( Color(1,0,0), 0.67*zmax );
	    palette.push_back( Color(0,0,0), zmax );
	}
	palette.norm();

	_colormap = new Colormap( range, histo2d->n(), histo2d->m(), histo2d->get_data() );
	_colormap->set_palette( palette );
	_colormap->set_interpolation( _interpolation );
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _colormap );

    } else if( _type == PARTICLE_DIAG_PLOT_HISTO1D ) {

	// Make XYGraph profile plot
	Histogram1D *histo1d = dynamic_cast<Histogram1D *>( _histo );
	if( !histo1d )
	    throw( Error( ERROR_LOCATION, "dynamic cast error" ) );
	if( _geom->geom_mode() == MODE_CYL && _diagx == DIAG_R ) {
	    _frame->set_axis_label( PLOT_AXIS_Y1, "Intensity (A/m^2)" );
	}
	std::vector<double> xdata;
	xdata.reserve( histo1d->n() );
	for( size_t a = 0; a < histo1d->n(); a++ )
	    xdata.push_back( histo1d->coord(a) );

	_profile = new XYGraph( xdata, histo1d->get_data() );
	_profile->set_color( Color(1,0,0) );
	_profile->set_line_style( XYGRAPH_LINE_SOLID );
	_profile->set_point_style( XYGRAPH_POINT_DISABLE );
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _profile );
    }

    if( (_type == PARTICLE_DIAG_PLOT_HISTO2D ||
	 _type == PARTICLE_DIAG_PLOT_SCATTER) &&
	((_diagx == DIAG_X && _diagy == DIAG_XP)  ||
	 (_diagx == DIAG_Y && _diagy == DIAG_YP)  ||
	 (_diagx == DIAG_R && _diagy == DIAG_RP)  ||
	 (_diagx == DIAG_Z && _diagy == DIAG_ZP)) && _ellipse_enable ) {

	// Plot emittance ellipse
	double a = _emit->rmajor();
	double b = _emit->rminor();
	double p = _emit->angle();
	double sinp = sin(p);
	double cosp = cos(p);
	std::vector<double> xd, yd;
	for( size_t i = 0; i < 100; i++ ) {
	    double t = 2.0*M_PI*i/99.0;
	    xd.push_back( _emit->xave()  + a*cos(t)*cosp - b*sin(t)*sinp );
	    yd.push_back( _emit->xpave() + a*cos(t)*sinp + b*sin(t)*cosp );
	}

	// Add ellipse xy graph
	_ellipse = new XYGraph( xd, yd );
	_ellipse->set_color( Color(0,0,0) );
	_ellipse->set_line_style( XYGRAPH_LINE_SOLID );
	_ellipse->set_point_style( XYGRAPH_POINT_DISABLE );
	_frame->add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, _ellipse );

	// Add emittance numbers to title
	std::stringstream ss;
	ss << "Emittance plot at " << coordinate_axis_string[_axis] << " = " << _level << " m\n"
	   << "\\alpha  = "   << _emit->alpha()   << ", "
	   << "\\beta  = "    << _emit->beta()    << " m/rad, "
	   << "\\gamma  = "   << _emit->gamma()   << " rad/m, "
	   << "\\epsilon  = " << _emit->epsilon() << " \\pi \\cdot m\\cdot rad";
	_frame->set_title( ss.str().c_str() );

    } else if( _type == PARTICLE_DIAG_PLOT_HISTO1D && 
	       (_diagx == DIAG_X || _diagx == DIAG_Y ||
		_diagx == DIAG_R || _diagx == DIAG_Z) ) {

	// Make title for profile
	std::stringstream ss;
	ss << "Profile plot at " << coordinate_axis_string[_axis] << " = " << _level << " m";
	_frame->set_title( ss.str().c_str() );

    } else {

	// Remove title
	_frame->set_title( "" );

    }
}







