#include <iostream>
#include <limits>
#include <algorithm>
#include "palette.hpp"



Palette::Entry::Entry( const Color &color, double val )
{
    _color = color;
    _val = val;
}


bool Palette::Entry::operator<( const Entry &e ) const
{
    return( _val < e._val );
}






Palette::Palette()
{
    _entries.push_back( Entry( Color( 1.0, 1.0, 1.0 ), 0.0 ) );
    _entries.push_back( Entry( Color( 0.0, 0.0, 0.0 ), 1.0 ) );
}


Palette::Palette( const std::vector<Entry> &entries )
{
    // Search minimum and maximum
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();
    for( size_t a = 0; a < entries.size(); a++ ) {
	if( entries[a]._val < min )
	    min = entries[a]._val;
	if( entries[a]._val > max )
	    max = entries[a]._val;
    }
    
    // Calculate offset and coef
    double offset = -min;
    double coef = 1.0/(max-min);
    if( max-min == 0.0 )
	coef = 1.0;

    // Clear old palette
    _entries.clear();

    // Add new entries
    _entries.reserve( entries.size() );
    for( size_t a = 0; a < entries.size(); a++ ) {
	_entries.push_back( Entry( entries[a]._color, coef*(offset+entries[a]._val) ) );
    }

    // Sort entries
    sort( _entries.begin(), _entries.end() );
}


void Palette::clear( void )
{
    // Clear old palette
    _entries.clear();
}


void Palette::push_back( const Color &color, double val )
{
    // Add new entry
    _entries.push_back( Entry( color, val ) );

    // Sort entries
    sort( _entries.begin(), _entries.end() );
}


void Palette::norm( void )
{
    // Search minimum and maximum
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();
    for( size_t a = 0; a < _entries.size(); a++ ) {
	if( _entries[a]._val < min )
	    min = _entries[a]._val;
	if( _entries[a]._val > max )
	    max = _entries[a]._val;
    }
    
    // Calculate offset and coef
    double offset = -min;
    double coef = 1.0/(max-min);
    if( max-min == 0.0 )
	coef = 1.0;

    // Renorm entries
    for( size_t a = 0; a < _entries.size(); a++ )
	_entries[a]._val = coef*(offset+_entries[a]._val);
}


Color Palette::operator()( double x ) const
{
    // If undefined
    if( _entries.size() == 0 )
	return( Color( 0, 0, 0 ) );
    else if( _entries.size() == 1 )
	return( _entries[0]._color );

    // If out of limits
    if( x < 0.0 )
	return( _entries[0]._color );
    else if( x > 1.0 )
	return( _entries[_entries.size()-1]._color );

    // Search correct index
    size_t a;
    for( a = 1; a < _entries.size(); a++ ) {
	if( x < _entries[a]._val )
	    break;
    }
    
    // Interpolate
    Color c = _entries[a-1]._color + 
	(x-_entries[a-1]._val)/(_entries[a]._val-_entries[a-1]._val)*
	(_entries[a]._color-_entries[a-1]._color);

    return( c );
}


void Palette::debug_print( void ) const
{
    std::cout << "**Palette\n";
    std::cout << "size = " << _entries.size() << "\n";
    for( size_t a = 0; a < _entries.size(); a++ ) {
	std::cout << "entries[" << a << "] = " 
		  << _entries[a]._val << " "
		  << _entries[a]._color[0] << " "
		  << _entries[a]._color[1] << " "
		  << _entries[a]._color[2] << " "
		  << _entries[a]._color[3] << "\n";
    }
}













