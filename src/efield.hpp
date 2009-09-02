/*! \file efield.hpp
 *  \brief Header file defining electric field class.
 */


#ifndef EFIELD_HPP
#define EFIELD_HPP 1


#include "vec3d.hpp"


/*! \brief Abstract base class for electric field.
 *
 *  Efield implementations can be made for 1) interpolating electric
 *  field mesh or 2) to use analytic expressions or 3) to interpolate
 *  and differentiate electric potential (EpotEfield does this). Of
 *  course other implementations are also possible.
 */
class Efield {

public:

    /*! \brief Virtual destructor.
     */
    virtual ~Efield() {}

    /*! \brief Operator for getting interpolated electric field value
     *  at \a x.
     */
    virtual Vec3D operator()( const Vec3D x ) const = 0;
};


#endif













