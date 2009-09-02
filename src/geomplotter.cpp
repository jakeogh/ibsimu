#include "geomplotter.hpp"


GeomPlotter::GeomPlotter( const Geometry *geom )
    : GeomPlot(&(Plotter::_frame), geom)
{
    // Don't use graph caches with non-interactive plotter
    disable_cache();
}


GeomPlotter::~GeomPlotter()
{

}


void GeomPlotter::build_plot( void )
{
    //ParticleDiagPlot::build_plot();
}
