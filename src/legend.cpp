#include "legend.hpp"


void LegendEntry::plot( cairo_t *cairo, double x, double y ) const
{

}


void LegendEntry::get_size( double &width, double &height ) const
{

}





void MultiEntryLegend::add_entry( const LegendEntry &entry )
{
    _entry.push_back( entry );
}


void MultiEntryLegend::plot( cairo_t *cairo, double x, double y ) const
{
    
}


void MultiEntryLegend::get_size( double &width, double &height ) const
{

}






void ColormapLegend::plot( cairo_t *cairo, double x, double y ) const
{
    
}


void ColormapLegend::get_size( double &width, double &height ) const
{

}


void ColormapLegend::set_height( double height )
{
    _height = height;
}














