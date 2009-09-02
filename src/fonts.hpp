/*! \file fonts.hpp
 *  \brief Header file font handling
 */


#ifndef FONTS_HPP
#define FONTS_HPP 1


#include <string>
#include <vector>
#include <cairo.h>
#include <fontconfig/fontconfig.h>
#include <ft2build.h>
#include FT_FREETYPE_H



/*! \brief %Font for %FontLib engine.
 */
class Font
{
    cairo_font_face_t   *_cff;
    FcPattern           *_pmatch;

    std::string          _family;
    cairo_font_slant_t   _slant;   /*!< \brief Slant of font.
				    *
				    * Either CAIRO_FONT_SLANT_NORMAL, 
				    * CAIRO_FONT_SLANT_ITALIC or CAIRO_FONT_SLANT_OBLIQUE. */
    cairo_font_weight_t  _weight;  /*!< \brief Weight of font.
				    *
				    * CAIRO_FONT_WEIGHT_NORMAL or CAIRO_FONT_WEIGHT_BOLD.
				    */

public:

    Font( FcPattern *pat );
    Font( const std::string &family );
    Font( const std::string &family, cairo_font_slant_t slant, 
	  cairo_font_weight_t weight );
    Font( const Font &font );

    Font &operator=( const Font &font );

    ~Font();

    std::string family( void ) const;
    cairo_font_slant_t slant( void ) const;
    cairo_font_weight_t weight( void ) const;
    cairo_font_face_t *font_face( void ) const;
    FcPattern *fcpattern( void ) const;

    void load( class FontLib &fontlib );

};


/*! \brief %Font engine using FreeType, FontConfig and cairographics.
 */
class FontLib
{

    FT_Library         _ft;
    FcConfig          *_fc;

    std::vector<Font>  _search;
    std::vector<Font>  _loaded;

    int get_font_matrix( cairo_t *cairo, cairo_matrix_t *matrix, cairo_matrix_t *orig_matrix );
    int process_substr( cairo_t *cairo, const std::string &str, cairo_text_extents_t *extents, 
			double x0, double y0, double &x, double &y );
    void process( cairo_t *cairo, const std::string &str, cairo_text_extents_t *extents, double &x, double &y );

public:

    FontLib();

    ~FontLib();

    /*! \brief Glyph symbol name entry.
     *
     *  Connects a user readable (LaTeX standard) name with a unicode code in utf-8.
     */
    struct Symbolname {
	const char *name;
	const char *ucode;
    };
    
    /*! \brief Combine extents.
     *
     *  Combine (extents1,x1,y1) and (extents2,x2,y2) to database
     * (extents1,x1,y1). The extents2 is taken to be latter and
     * therefore advance is defined to be from origo of the first data
     * to the advance of the second data. */
    static void combine_extents( cairo_text_extents_t *extents1, double x1, double y1,
				 const cairo_text_extents_t *extents2, double x2, double y2 );

    /*! \brief Chart of glyph symbol names.
     */
    static const Symbolname symbols[];


    FcConfig *fc( void ) { return( _fc ); }


    void push_auto_search_font( const std::string &family );
    int pop_auto_search_font( void );


    
    std::string family( void ) const;
    cairo_font_slant_t slant( void ) const;
    cairo_font_weight_t weight( void ) const;
    cairo_font_face_t *font_face( void ) const;
    FcPattern *fcpattern( void ) const;


    void push_font( FcPattern *pat );
    void push_font( const std::string &family, cairo_font_slant_t slant, 
		    cairo_font_weight_t weight );
    int pop_font( void );


    void text_extents( cairo_t *cairo, const std::string &str, cairo_text_extents_t *extents );

    /*! \brief Draw piece of text at (x,y)
     *
     *  The (x,y) are updated according to cursor advance.
     */
    void draw_text( cairo_t *cairo, const std::string &str, double &x, double &y );

};


extern FontLib fontlib;


#endif













