/*! \file hardcopy.hpp
 *  \brief Header file for simulation hardcopy.
 */


#ifndef GTK_HARDCOPY_HPP
#define GTK_HARDCOPY_HPP 1


#include <gtk/gtk.h>
#include <png.h>
#include "frame.hpp"




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













