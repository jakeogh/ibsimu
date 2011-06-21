/*! \file legend.cpp
 *  \brief Plot legends
 */

/* Copyright (c) 2005-2009,2011 Taneli Kalvas. All rights reserved.
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

#include "legend.hpp"


// Legend sizes are relative to font size
#define LEGEND_SAMPLE_WIDTH 3.0
#define LEGEND_SAMPLE_HEIGHT 1.0
#define LEGEND_SAMPLE_LINESKIP 1.5
#define LEGEND_SAMPLE_SEPARATION 0.3
#define LEGEND_MARGIN 0.5


void LegendEntry::plot( cairo_t *cairo, double x, double y )
{
    double scale = _label.get_font_size();
    _graph.plot_sample( cairo, x, y, scale*LEGEND_SAMPLE_WIDTH, scale*LEGEND_SAMPLE_HEIGHT );
    _label.set_location( x + scale*(LEGEND_SAMPLE_WIDTH + LEGEND_SAMPLE_SEPARATION), 
			 y );
    _label.draw( cairo );
}


void LegendEntry::get_size( cairo_t *cairo, double &width, double &height ) const
{
    double bbox[4];
    double scale = _label.get_font_size();
    _label.get_bbox( cairo, bbox );
    width = scale*(LEGEND_SAMPLE_WIDTH + LEGEND_SAMPLE_SEPARATION)
	+ (bbox[2]-bbox[0]);
    height = scale*LEGEND_SAMPLE_HEIGHT;
}


void LegendEntry::set_font_size( double fontsize )
{
    _label.set_font_size( fontsize );
}


/* ********************************************************************* */



MultiEntryLegend::MultiEntryLegend()
    : _fontsize(12.0)
{

}


void MultiEntryLegend::add_entry( LegendEntry *entry )
{
    if( entry )
	entry->set_font_size( _fontsize );
    _entry.push_back( entry );
}


void MultiEntryLegend::clear_entries( void )
{
    _entry.clear();
}


void MultiEntryLegend::plot( cairo_t *cairo, double x, double y )
{
    uint32_t k = 0;
    for( int32_t i = _entry.size()-1; i >= 0; i-- ) {
	if( _entry[i] ) {
	    _entry[i]->plot( cairo, 
			     x + _fontsize*LEGEND_MARGIN, 
			     y - _fontsize*(LEGEND_SAMPLE_LINESKIP*k+LEGEND_MARGIN) );
	    k++;
	}
    }
}


void MultiEntryLegend::get_size( cairo_t *cairo, double &width, double &height ) const
{
    width = height = 0.0;
    for( uint32_t i = 0; i < _entry.size(); i++ ) {
	double w, h;
	if( _entry[i] ) {
	    _entry[i]->get_size( cairo, w, h );
	    if( w > width )
		width = w;
	    if( i == _entry.size()-1 )
		height += _fontsize*LEGEND_SAMPLE_HEIGHT;
	    else
		height += _fontsize*LEGEND_SAMPLE_LINESKIP;
	}
    }

    height += _fontsize*2.0*LEGEND_MARGIN;
    width += _fontsize*2.0*LEGEND_MARGIN;
}


void MultiEntryLegend::set_font_size( double fontsize )
{
    _fontsize = fontsize;
    for( uint32_t i = 0; i < _entry.size(); i++ ) {
	if( _entry[i] )
	    _entry[i]->set_font_size( fontsize );
    }
}


/* ********************************************************************* */



void ColormapLegend::plot( cairo_t *cairo, double x, double y )
{
    
}


void ColormapLegend::get_size( cairo_t *cairo, double &width, double &height ) const
{

}


void ColormapLegend::set_height( double height )
{
    _height = height;
}




















