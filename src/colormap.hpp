/*! \file colormap.hpp
 *  \brief Header file for plot colormap
 */


#ifndef COLORMAP_HPP
#define COLORMAP_HPP 1


#include <cairo.h>
#include <vector>
#include "palette.hpp"
#include "graph.hpp"
#include "coordmapper.hpp"



enum interpolation_e {
    INTERPOLATION_CLOSEST = 0,
    INTERPOLATION_BILINEAR,
    INTERPOLATION_BICUBIC
};


enum zscale_e {
    ZSCALE_LINEAR = 0,
    ZSCALE_LOG,
    ZSCALE_RELLOG
};


/*! \brief Class for colormap type plots.
 *
 *  Implementation of %Graph.
 */
class Colormap : public Graph {

    Palette                _palette;       /*!< \brief Palette for plotting. */

    interpolation_e        _interpolation; /*!< \brief Interpolation mode. */
    zscale_e               _zscale;        /*!< \brief zscale mode. */
    
    double                 _zmin;          /*!< \brief Minimum zval. */
    double                 _zmax;          /*!< \brief Maximum zval. */

    double                 _datarange[4];  /*!< \brief Data ranges: xmin, ymin, xmax, ymax. */
    size_t                 _n;             /*!< \brief Size of data-array in x-direction. */
    size_t                 _m;             /*!< \brief Size of data-array in y-direction. */

    std::vector<double>    _f;             /*!< \brief Functuion value data, y major order. */

public:

    /*! \brief Default constructor for empty colormap graph.
     */
    Colormap();

    /*! \brief Constructor for basic graph with defined data.
     *
     *  Data is defined as \a n by \a m array of data, where x and y
     *  ranges are defined in datarange in order xmin, ymin, xmax,
     *  ymax. Z-values are defined in vector \a data in y major
     *  order. Internal copies of the data from data is made.
     */
    Colormap( const double datarange[4], size_t n, size_t m, 
	      const std::vector<double> &data );

    /*! \brief Destructor.
     */
    virtual ~Colormap() {}

    /*! \brief Set interpolation mode.
     */
    void set_interpolation( interpolation_e interpolation );

    /*! \brief Set zscale mode.
     */
    void set_zscale( zscale_e zscale );

    /*! \brief Plot graph with cairo.
     *
     *  Plot the graph using \a cairo and coordinate mapper \a
     *  cm. The visible range of plot is given in array \a range in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] );

    /*! \brief Get bounding box of graph.
     *
     *  Returns the bounding box of the graph in array \a bbox in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void get_bbox( double bbox[4] );

    /*! \brief Set colormap palette.
     */
    void set_palette( const Palette &palette );

    /*! \brief Get zrange for colormap plot.
     */
    void get_zrange( double &min, double &max ) const;
};


#endif















