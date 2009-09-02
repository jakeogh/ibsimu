#include "fielddiagplotter.hpp"


FieldDiagPlotter::FieldDiagPlotter( const Geometry *geom )
    : FieldDiagPlot(&(Plotter::_frame),geom)
{

}


FieldDiagPlotter::~FieldDiagPlotter()
{
}


void FieldDiagPlotter::build_plot( void )
{
    FieldDiagPlot::build_plot();
}


