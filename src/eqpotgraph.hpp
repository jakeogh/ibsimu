/*! \file eqpotgraph.hpp
 *  \brief Header file for equipotential line plotter.
 */


#ifndef EQPOTPLOT_HPP
#define EQPOTPLOT_HPP 1


#include <vector>
#include "geometry.hpp"
#include "scalarfield.hpp"
#include "graph3d.hpp"
#include "color.hpp"


/*! \brief Equipotential line plot.
 *
 *  Class for constructing and drawing equipotential line plots.
 */
class EqPotGraph : public Graph3D {

    struct Line {
	double x[4];              /* Coordinates of line segment (x1,y1,x2,y2) */
	
	Line( double _x1, double _y1, double _x2, double _y2 ) { 
	    x[0] = _x1; 
	    x[1] = _y1;
	    x[2] = _x2; 
	    x[3] = _y2;
	}

	double &operator[]( int i ) { return( x[i] ); }
	const double &operator[]( int i ) const { return( x[i] ); }
    };

    struct EqPotLines {
	double            pot;    /* Potential value for equipotential line. */
	std::vector<Line> x;      /* Line segment coordinates. */
	
	EqPotLines( double pot ) : pot(pot) {}
    };

    Color                               _color;
    const ScalarField                   _epot;
    const Geometry                     &_g;
    bool                                _data_built;

    std::vector<double>                 _eqlines_manual;
    size_t                              _eqlines_auto;
    std::vector<EqPotLines *>           _lines;

    view_e                              _oview;
    double                              _olevel;

    bool                                _cache;


    bool eqline_exists( double pot1, signed char sol1, 
			double pot2, signed char sol2, 
			double pot ) const;
    void build_data( void );

public:

    /*! \brief Constructor for equipotential line plot.
     *
     *  Makes a plot object for plotting equipotential data from
     *  scalarfield \a field in geometry \a g.
     */
    EqPotGraph( const ScalarField &epot, const Geometry &g );

    /*! \brief Destructor,
     */
    virtual ~EqPotGraph();

    /*! \brief Disable internal cache.
     *
     *  Makes equipotential lines to be calculated at every plot().
     */
    void disable_cache( void );

    /*! \brief Add manual equipotential lines to be plotted at
     *  specified potentials.
     */
    void set_eqlines_manual( const std::vector<double> &pot );

    /*! \brief Set \a N automatic equipotential lines to be plotted
     *  between minimum potential and maximum potentials.
     */
    void set_eqlines_auto( size_t N );

    /*! \brief Plot drawable with cairo.
     *
     *  Plot the drawable using \a cairo and coordinate mapper \a
     *  cm. The visible range of plot is given in array \a range in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] );

    /*! \brief Get bounding box of drawable.
     *
     *  Returns the bounding box of the drawable in array \a bbox in
     *  order xmin, ymin, xmax, ymax.
     */
    virtual void get_bbox( double bbox[4] );
};


#endif













