/*! \file solid.hpp
 *  \brief Header file defining class solid.
 */


#ifndef SOLID_HPP
#define SOLID_HPP 1


#include <iostream>
#include "vec3d.hpp"


#define FILEID_FUNCSOLID 2001
#define FILEID_CSGSOLID 2002


/*! \brief Abstract solid class.
 *
 *  %Solid class holds the definition for one solid. %Solid class is
 *  an abstract base class. Different implementation may exist.
 */
class Solid {
    
public:

    /*! \brief Virtual destructor.
     */
    virtual ~Solid() {}

    /*! \brief Return if point x is inside solid.
     */
    virtual bool inside( const Vec3D &x ) const = 0;

    /*! \brief Prints internal data to std::cout.
     */
    virtual void debug_print( void ) const {}

    /*! \brief Saves solid data to stream.
     */
    virtual void save( std::ostream &s ) const = 0;
};


#endif













