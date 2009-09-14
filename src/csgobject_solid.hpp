/*! \file csgobject_solid.hpp
 *  \brief Header file defining class CSGObjectSolid.
 */


#ifndef CSGOBJECT_SOLID_HPP
#define CSGOBJECT_SOLID_HPP 1


#include <iostream>
#include <csg.hpp>
#include "solid.hpp"


/*! \brief CSGObject solid class.
 *
 *  %CSGObjectSolid is an implementation of Solid using constructive
 *  solid geometry objects.
 */
class CSGObjectSolid : public Solid {

    CSGObject *_object;

public:

    /*! \brief Constructor for Solid using %CSGObject \a object for
     *  defining the geometry.
     *
     *  A pointer to the %CSGObject is saved in the
     *  %CSGObjectSolid. The user must take care of calling the
     *  destructor of the %CSGObject after use.
     */
    CSGObjectSolid( CSGObject *object ) : _object(object) {}

    /*! \brief Constructor for loading solid data from a file.
     */
    CSGObjectSolid( std::istream &s );

    /*! \brief Destructor.
     */
    ~CSGObjectSolid() {}

    /*! \brief Return if point x is inside funcsolid.
     */
    bool inside( const Vec3D &x ) const;

    /*! \brief Prints internal data to std::cout.
     */
    void debug_print( void ) const;

    /*! \brief Saves solid data to stream.
     */
    void save( std::ostream &fout ) const;
};


#endif















