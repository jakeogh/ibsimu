/*! \file xygraph.hpp
 *  \brief Header file for plot XYGraph
 */


#ifndef XYGRAPH_HPP
#define XYGRAPH_HPP 1


#include <cairo.h>
#include <vector>
#include "color.hpp"
#include "graph.hpp"
#include "coordmapper.hpp"


enum line_style_e {
    XYGRAPH_LINE_DISABLE = 0,
    XYGRAPH_LINE_SOLID
};


enum point_style_e {
    XYGRAPH_POINT_DISABLE = 0,
    XYGRAPH_POINT_CIRCLE
};


/*! \brief Class for XY-type simple graph plots. 
 *
 *  Implementation of %Graph.
 */
class XYGraph : public Graph {

    double                 _linewidth;
    Color                  _color; 
    line_style_e           _linestyle;
    point_style_e          _pointstyle;
    bool                   _point_filled;
    double                 _point_scale;
    
    std::vector<double>    _xdata;
    std::vector<double>    _ydata;

public:

    /*! \brief Default constructor for empty graph.
     */
    XYGraph();

    /*! \brief Constructor for basic graph with datapoints \a xdata
     *  and \a ydata.
     *
     *  Internal copies of the data from xdata and ydata are made.
     */
    XYGraph( const std::vector<double> &xdata, 
	     const std::vector<double> &ydata );

    /*! \brief Destructor.
     */
    virtual ~XYGraph() {}

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

    /*! \brief Set new data arrays.
     */
    void set_data( const std::vector<double> &xdata, 
		   const std::vector<double> &ydata );

    /*! \brief Set line width.
     */
    void set_line_width( double linewidth );

    /*! \brief Set graph color.
     */
    void set_color( const Color &color );

    /*! \brief Set line style.
     */
    void set_line_style( line_style_e linestyle );

    /*! \brief Set point style.
     */
    void set_point_style( point_style_e pointstyle, bool filled = true, double scale = 1.0 );
};


#endif













