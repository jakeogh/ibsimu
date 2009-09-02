/*! \file func_solid.hpp
 *  \brief Header file defining class FuncSolid.
 */


#ifndef FUNC_SOLID_HPP
#define FUNC_SOLID_HPP 1


#include <iostream>
#include "solid.hpp"


/*! \brief Function solid class.
 *
 *  FuncSolid class holds the definition for one solid defining
 *  C-function. This solid implementation suffers from the inability
 *  of saving to file. If FuncSolid is constructed from stream the
 *  function pointer inside is set to NULL and error is thrown if
 *  function is evaluated using inside().
 */
class FuncSolid : public Solid {

    bool (*_func)(double,double,double);

public:

    /*! \brief Constructor.
     */
    FuncSolid( bool (*func)(double,double,double) ) : _func(func) {}

    /*! \brief Constructor for loading solid data from a file.
     */
    FuncSolid( std::istream &s );

    /*! \brief Destructor.
     */
    ~FuncSolid() {}

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













