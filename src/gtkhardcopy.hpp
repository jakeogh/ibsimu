/*! \file gtkhardcopy.hpp
 *  \brief Dialog window for producing hard copies
 */

/* Copyright (c) 2005-2009 Taneli Kalvas. All rights reserved.
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

#ifndef GTKHARDCOPY_HPP
#define GTKHARDCOPY_HPP 1


#include <gtk/gtk.h>
#include <png.h>
#include "frame.hpp"




/*! \brief Interactive dialog for producing hardcopies.
 */
class GTKHardcopy
{
    GtkWidget        *_window;
    Frame            *_frame;

    double            _aspect;
    size_t            _width;
    size_t            _height;

    GtkWidget        *_spinx;
    GtkWidget        *_spiny;

    static void spinx_signal( GtkSpinButton *spinbutton,
			      gpointer object );
    static void spiny_signal( GtkSpinButton *spinbutton,
			      gpointer object );
    static int type_from_extension( const char *filename );
    static void ensure_extension( std::string &filename, 
				  const std::string &extension );
    static void treeview_changed( GtkTreeSelection *selection,
				  gpointer userdata );
    
    void spinx( void );
    void spiny( void );

    void get_image_size( cairo_surface_t *p_surface, 
			 int &width, int &height );
    static void unpremultiply_data( png_structp png, 
				    png_row_infop row_info, 
				    png_bytep data );
    void write_to_png( cairo_surface_t *p_surface, 
		       int width, int height, 
		       const char *filename );
    void write_png( const char *filename );


    void write_eps( const char *filename );


    void write_svg( const char *filename );


    void write_pdf( const char *filename );
public:

    GTKHardcopy( GtkWidget *window, Frame *frame, size_t width, size_t height );

    ~GTKHardcopy();

    void run( void );

};


#endif


















