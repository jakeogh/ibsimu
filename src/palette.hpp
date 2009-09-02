/*! \file palette.hpp
 *  \brief Header file for palette.
 */


#ifndef PALETTE_HPP
#define PALETTE_HPP 1


#include <vector>
#include "color.hpp"


/*! \brief Class for palette definition.
 *
 *  Palette is an object that contains a list of colors and
 *  corresponding values. The colors are interpolated linearly between
 *  the defined points for a smooth color palette. The palette values
 *  are normed so that the end points of palette have values 0.0 and
 *  1.0.
 */
class Palette {

public:

    struct Entry {
	Color  _color; /*!< \brief Entry colors */
	double _val;   /*!< \brief Entry values. */

	Entry( const Color &color, double val );

	bool operator<( const Entry &e ) const;
    };

private:

    std::vector<Entry> _entries;    /*!< \brief Palette entries. */

public:

    /*! \brief Default constructor for default black and white palette.
     */
    Palette();

    /*! \brief Constructor for defined palette.
     *
     *  The values are normed to range from 0.0 to 1.0.
     */
    Palette( const std::vector<Entry> &entries );

    /*! \brief Return the interpolated color value from palette.
     *
     *  Makes and interpolated color value at value \a x. If palette
     *  has no colors, black will be returned. If palette has one
     *  color, that color will be returned. With two or more colors
     *  the interpolation can be done correctly.
     */
    Color operator()( double x ) const;

    /*! \brief Clear current palette.
     *
     *  Leaves palette with no colors. Used with push_back() to build
     *  new palettes on-line.
     */
    void clear( void );

    /*! \brief Pushes new entry to palette.
     *
     *  Palette entries are automatically sorted. The palette won't be
     *  normed in range. This has to be manually done by calling
     *  norm() after adding palette entries.
     */
    void push_back( const Color &color, double val );

    /*! \brief Normalize palette entries.
     *
     *  Normalize palette to range from 0.0 to 1.0.
     */
    void norm( void );

    void debug_print( void ) const;
};


#endif













