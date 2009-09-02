/*! \file ruler.hpp
 *  \brief Header file for plot ruler class
 */


#ifndef RULER_HPP
#define RULER_HPP 1


#include <cairo.h>
#include <vector>
#include <string>
#include <iostream>
#include "color.hpp"
#include "label.hpp"
#include "coordmapper.hpp"


/*! \brief Class for coordinate axis ruler.
 *
 *  Ruler is a class for drawing coordinate axis for frame or palette
 *  colormap legend for example.
 */
class Ruler {

    /*! \brief One coordinate tic with label and location.
     */
    struct Tic {
	double               _x;            /*!< \brief Tic coordinate. */
	Label                _label;        /*!< \brief Tic label. */
	
	/*! \brief Constructor for coordinate tic.
	 */
	Tic( double x, const std::string &label ) : _x(x), _label(label) {}
    };

    Color                    _color;        /*!< \brief Ruler color. */
    double                   _ticlen_in;    /*!< \brief Tic length in. */
    double                   _ticlen_out;   /*!< \brief Tic length out. */
    double                   _labelspace;   /*!< \brief Space from tic to ticlabel. */

    double                   _range[2];     /*!< \brief Coordinate range. */
    bool                     _autorange[2]; /*!< \brief Autorange enable. */

    double                   _endpt[4];     /*!< \brief Ruler endpoints (x1,y1,x2,y2). */

    double                   _fontsize;     /*!< \brief Ruler font size. */
    Label                    _axislabel;    /*!< \brief Ruler axislabel. */

    bool                     _label_enabled;
    bool                     _indir;        /*!< \brief In direction. True for x1 and y2 axes. */
    std::vector<Tic>         _tic;          /*!< \brief Tics. */
    
    int                      _cind;         /*!< \brief Coordinate index, 0 or 1. */


    static bool tic_labels_bbox_crash_x( const double bbox[4], const double obbox[4] );
    static bool tic_labels_bbox_crash_y( const double bbox[4], const double obbox[4] );
    bool add_tic( double x, cairo_t *cairo, const Coordmapper1D &cm, 
		  bool ruler_tic_bbox_test, double &maxsize, double obbox[4] );

public:

    /*! \brief Default contructor.
     */
    Ruler();

    /*! \brief Contructor for defined coordinate index.
     */
    Ruler( int cind );

    /*! \brief Copy contructor.
     */
    Ruler( const Ruler &ruler );

    /*! \brief Assignment operator.
     */
    Ruler &operator=( const Ruler &ruler );

    /*! \brief Destructor.
     */
    ~Ruler() {}

    /*! \brief Copy tics and from other ruler.
     */
    void copy_tics( const Ruler &ruler );

    /*! \brief Set font size.
     */
    void set_font_size( double size );

    /*! \brief Set ruler color.
     */
    void set_color( const Color &color );

    /*! \brief Set tic lengths
     */
    void set_ticlen( double inlen, double outlen );

    /*! \brief Get ruler autoranging status
     */
    void set_autorange( bool autorange_min, bool autorange_max );

    /*! \brief Get ruler autoranging status
     */
    void get_autorange( bool &autorange_min, bool &autorange_max ) const;

    /*! \brief Set ruler ranges.
     *
     *  Set ruler ranges to finite values. If only one range is either
     *  infinite or NaN, it is replaced by a number to get a total
     *  range span of 1. If both numbers are infinite or NaN, the
     *  minimum will be set to 0.0 and maximum to 1.0.
     */
    void set_ranges( double min, double max );

    /*! \brief Get current ruler ranges.
     *
     *  Returns finite values used for ruler ranges.
     */
    void get_ranges( double &min, double &max ) const;

    /*! \brief End points in pixels for ruler.
     */
    void set_endpoints( double x1, double y1, double x2, double y2 );

    /*! \brief Set axis label.
     */
    void set_axis_label( const std::string &label );

    /*! \brief Enable axis title and tic labels.
     */
    void enable_labels( bool enable );

    /*! \brief Set coordinate index.
     *
     *  Set \a cind to 0 for x and 1 for y.
     */
    void set_coord_index( int cind );

    /*! \brief Set in direction for frame ruler.
     */
    void set_indir( bool ccw );

    /*! \brief Calculate ruler ranges, set tics and coordmapper.
     *
     *  If ruler_tic_bbox_test is set to true, the bounding boxes of
     *  the tic labels are tested against overlapping and tics are
     *  rebuilt with larger tic spacing to prevent this.
     */
    void calculate( cairo_t *cairo, Coordmapper1D &cm, bool ruler_tic_bbox_test );

    /*! \brief Draw ruler with cairo.
     *
     *  If \a recalculate is false, the ruler will be drawn with it's
     *  current settings. This is for frame to prevent excess
     *  recalculation of rulers.
     */
    void draw( cairo_t *cairo, Coordmapper1D &cm, bool recalculate = true );

    /*! \brief Get bounding box of ruler.
     *
     *  If \a recalculate is false, the ruler bounding box will
     *  correspond to it's current settings. This is for frame to
     *  prevent excess recalculation of rulers.
     */
    void get_bbox( cairo_t *cairo, double bbox[4], Coordmapper1D &cm, bool recalculate = true );

    /*! \brief Debug print to stream.
     */ 
    void debug_print( std::ostream &os ) const ;
};


#endif













