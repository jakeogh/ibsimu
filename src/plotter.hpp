#ifndef PLOTTER_HPP
#define PLOTTER_HPP 1


#include <png.h>
#include <cairo.h>
#ifdef CAIRO_HAS_SVG_SURFACE
#include <cairo-svg.h>
#endif
#if CAIRO_HAS_PS_SURFACE
#include <cairo-ps.h>
#endif
#ifdef CAIRO_HAS_PDF_SURFACE
#include <cairo-pdf.h>
#endif
#include "frame.hpp"


/*! \brief Non-interactive plotter.
 *
 *  Plotter skeleton for building plots non-interactively.
 */
class Plotter {

    size_t _width;
    size_t _height;

#ifdef CAIRO_HAS_PNG_FUNCTIONS
    static void png_get_image_size( cairo_surface_t *p_surface, 
				    int &width, int &height );
    static void png_unpremultiply_data( png_structp png, 
					png_row_infop row_info, 
					png_bytep data );
    static void write_to_png( cairo_surface_t *p_surface, 
			      int width, int height, 
			      const char *filename );
#endif

    virtual void build_plot( void ) = 0;

protected:

    Frame _frame;

    /*! \brief Constructor for plotter.
     *
     *  Not intended to be used on its own.
     */
    Plotter();

    /*! \brief Destructor for plotter.
     */
    virtual ~Plotter();

public:

    /*! \brief Set size of plot.
     */
    void set_size( size_t width, size_t height ) {
	_width = width;
	_height = height;
    } 

    /*! \brief Set font size for plot.
     */
    void set_font_size( size_t size );

    /*! \brief Set ranges of plot in x- and y-directions.
     */
    void set_ranges( double xmin, double ymin, double xmax, double ymax );

#ifdef CAIRO_HAS_PNG_FUNCTIONS
    /*! \brief Make a plot to a PNG-file.
     */
    void plot_png( const std::string &filename );
#endif

#ifdef CAIRO_HAS_PS_SURFACE
    /*! \brief Make a plot to a EPS-file.
     */
    void plot_eps( const std::string &filename );
#endif

#ifdef CAIRO_HAS_PDF_SURFACE
    /*! \brief Make a plot to a PDF-file.
     */
    void plot_pdf( const std::string &filename );
#endif

#ifdef CAIRO_HAS_SVG_SURFACE
    /*! \brief Make a plot to a SVG-file.
     */
    void plot_svg( const std::string &filename );
#endif

};



#endif
