#include "particlediagplotter.hpp"


ParticleDiagPlotter::ParticleDiagPlotter( const Geometry *geom, const ParticleDataBase *pdb, 
					  coordinate_axis_e axis, double level, 
					  particle_diag_plot_type_e type,
					  trajectory_diagnostic_e diagx, trajectory_diagnostic_e diagy )
    : ParticleDiagPlot(&(Plotter::_frame), geom, pdb, axis, level, type, diagx, diagy)
{

}


ParticleDiagPlotter::~ParticleDiagPlotter()
{

}


void ParticleDiagPlotter::build_plot( void )
{
    ParticleDiagPlot::build_plot();
}
