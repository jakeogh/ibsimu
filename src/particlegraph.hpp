/*! \file particlegraph.hpp
 *  \brief Header file for particle trajectory plotter.
 */


#ifndef PARTICLEPLOT_HPP
#define PARTICLEPLOT_HPP 1


#include <cairo.h>
#include <vector>
#include "color.hpp"
#include "graph3d.hpp"
#include "coordmapper.hpp"
#include "geometry.hpp"
#include "particledatabase.hpp"
#include "lineclip.hpp"


/*! \brief Class for drawing particle trajectories.
 *
 *  Implementation of %Graph3D.
 */
class ParticleGraph : public Graph3D {
    
    const Geometry         &_g;          /*!< \brief Reference to simulation geometry. */
    const ParticleDataBase &_pdb;        /*!< \brief Reference to particle database. */
    int                     _particlediv;/*!< \brief Particle plot divisor. */

    std::vector<Color>      _color;      /*!< \brief Colors for trajectories. */

    double                  _ox[5];      /*!< \brief Workspace for particleplot_draw_cruve() */
    size_t                  _coordsize;  /*!< \brief Size of array _coord divided by two */
    double                 *_coord;      /*!< \brief Workspace for particleplot_draw_curve() */
    bool                    _qm_discr;   /*!< \brief q/m discriminator enable, default true */

    void get_point( const Coordmapper *cm, double *coord, double s, 
		    double Ax, double Bx, double Cx, double Dx, 
		    double Ay, double By, double Cy, double Dy ) const;

    void draw_linear( const Coordmapper *cm, LineClip &lc, 
		      double x[5], bool first ) const;

    void draw_curve( const Coordmapper *cm, LineClip &lc, 
		     double x[5], bool first );

public:

    /*! \brief Constructor for particle plotter.
     */
    ParticleGraph( const Geometry &g, const ParticleDataBase &pdb,
		  int particlediv = 10, bool qm_discr = true );

    /*! \brief Destructor.
     */
    virtual ~ParticleGraph();

    /*! \brief Set particle divisor.
     *
     *  Zero for no plotting, one for plotting every particle, two for
     *  plotting every second particle, three for plotting every
     *  third particle, etc. Defaults to 10.
     */
    void set_particle_div( size_t particle_div );

    /*! \brief Enable q/m discretation
     *
     *  If q/m discretation is enabled, different q/m particles will
     *  be plotted with different colors. Default is enabled.
     */
    void set_qm_discretation( bool qm_discr );

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

    /*! \brief Add a color to the list of trajectory colors.
     */
    void add_color( const Color &color );

    /*! \brief Clear the list of trajectory colors.
     */
    void clear_colors( void );

};



#endif













