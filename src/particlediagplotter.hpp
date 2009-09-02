#include "plotter.hpp"
#include "particlediagplot.hpp"


/*! \brief Non-interactive particle diagnostic plotter.
 *
 *  This class ties together Plotter, which provides basic graphics
 *  features including frame control, and ParticleDiagPlot, which
 *  makes the particle diagnostic plot.
 */
class ParticleDiagPlotter : public Plotter, public ParticleDiagPlot {

    virtual void build_plot( void );

public:

    /*! \brief Constructor for particle diagnostic plotter.
     */
    ParticleDiagPlotter( const Geometry *geom, const ParticleDataBase *pdb, 
			 coordinate_axis_e axis, double level, 
			 particle_diag_plot_type_e type,
			 trajectory_diagnostic_e diagx, trajectory_diagnostic_e diagy );

    /*! \brief Destructor for particle diagnostic plotter.
     */
    ~ParticleDiagPlotter();

};
