/*! \file legend.hpp
 *  \brief Header file for legend.
 */


#ifndef LEGEND_HPP
#define LEGEND_HPP 1


#include <vector>
#include <string>
#include "graph.hpp"
#include "colormap.hpp"


/*! \brief Class for legend entry.
 *
 *  The legend entries contain a reference to the drawable so that if
 *  the style used in the drawable is changed, the sample in legend is
 *  automatically changed.
 */
class LegendEntry {
    
    Graph     &_drawable;
    std::string   _label;

public:
    
    /*! \brief Contructor for legend entry.
     */
    LegendEntry( Graph &drawable, std::string &label ) 
	: _drawable(drawable), _label(label) {}

    /*! \brief Copy constructor.
     */
    LegendEntry( const LegendEntry &le ) 
	: _drawable(le._drawable), _label(le._label) {}

    /*! \brief Destructor.
     */
    ~LegendEntry() {}

    /*! \brief Assignment operator.
     */
    LegendEntry &operator=( const LegendEntry &le ) {
	_drawable = le._drawable;
	_label = le._label;
	return( *this );
    }

    /*! \brief Plot legend entry.
     */
    void plot( cairo_t *cairo, double x, double y ) const;

    /*! \brief Get size of legend entry.
     */
    void get_size( double &width, double &height ) const;

};


/*! \brief Base class for legend definition.
 *
 *  Legend is an object that contains a key to the plot styles used in
 *  Graphs. The key contains a sample of the plot style used and a
 *  corresponding text label. 
 *
 *  The Colormap legend is a special case because in addition to the
 *  plot style, the plot z-range is shown in the legend.
 */
class Legend {

public:

    /*! \brief Default constructor for legend.
     */
    Legend() {}

    /*! \brief Virtual destructor.
     */
    virtual ~Legend() {}    

    /*! \brief Plot legend at (x,y).
     */
    virtual void plot( cairo_t *cairo, double x, double y ) const = 0;

    /*! \brief Get size of legend.
     */
    virtual void get_size( double &width, double &height ) const = 0;
};


/*! \brief Legend for presenting plot styles.
 */
class MultiEntryLegend {

    std::vector<LegendEntry> _entry;    /*!< \brief Legend entries. */

public:

    /*! \brief Default constructor for legend.
     */
    MultiEntryLegend() {}
    
    /*! \brief Virtual destructor.
     */
    virtual ~MultiEntryLegend() {}    

    /*! \brief Plot legend at (x,y).
     */
    virtual void plot( cairo_t *cairo, double x, double y ) const;

    /*! \brief Get size of legend.
     */
    virtual void get_size( double &width, double &height ) const;

    /*! \brief Add entry to legend.
     */
    void add_entry( const LegendEntry &entry );
};


/*! \brief Legend for presenting colormap key.
 */
class ColormapLegend {

    double      _height;
    Colormap   &_colormap;


public:

    /*! \brief Default constructor for legend.
     */
    ColormapLegend( Colormap &colormap ) : _height(0.0), _colormap(colormap) {}

    /*! \brief Virtual destructor.
     */
    virtual ~ColormapLegend() {}

    /*! \brief Plot legend at (x,y).
     */
    virtual void plot( cairo_t *cairo, double x, double y ) const;

    /*! \brief Get size of legend.
     */
    virtual void get_size( double &width, double &height ) const;

    /*! \brief Set height of legend.
     */
    void set_height( double height );
};


#endif













