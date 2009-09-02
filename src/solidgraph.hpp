/*! \file solidgraph.hpp
 *  \brief Header file for solid drawable.
 */


#ifndef SOLIDPLOT_HPP
#define SOLIDPLOT_HPP 1


#include <vector>
#include "geometry.hpp"
#include "graph3d.hpp"
#include "color.hpp"


/*! \brief A 2D cut view of the geometry solids.
 *
 *  Class for constructing and plotting a view of the geometry
 *  solids. The view data is stored inside the object in a cache to
 *  speed up more frequent use (in interactive plotter).
 */
class SolidGraph : public Graph3D {

    struct Point {
	double x[2];

	Point( double _x, double _y ) { x[0] = _x; x[1] = _y; }
	
	double &operator[]( int i ) { return( x[i] ); }
	const double &operator[]( int i ) const { return( x[i] ); }
    };

    struct SolidPoints {
	int                N;    /* Solid number of electrode */
	std::vector<Point> p;    /* Coordinate points of electrode boundary */
	
	SolidPoints( int N ) : N(N) {}
    };

    Color                                _color;
    const Geometry                      &_g;
    std::vector<SolidPoints *>           _solid;

    view_e                               _oview;
    double                               _olevel;
    
    bool                                 _cache;

    bool is_edge( int N, const int i[3] ) const;
    int get_mesh( const int i[3], int offsetx, int offsety ) const;
    void build_solid( SolidPoints *solid, const int i[3], char *done, int last, int N );
    void build_data( void );
    void clear_data( void );
    
public:

    /*! \brief Constructor for %SolidGraph drawable from geometry \a g.
     */
    SolidGraph( const Geometry &g );

    /*! \brief Destructor.
     */
    virtual ~SolidGraph();

    /*! \brief Disable internal cache.
     *
     *  Makes solid boundaries to be calculated at every plot().
     */
    void disable_cache( void );

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













