/*! \file solver.hpp
 *  \brief Header file defining class solver.
 */


#ifndef SOLVER_HPP
#define SOLVER_HPP 1


#include <iostream>
#include "geometry.hpp"
#include "problem.hpp"
#include "matrix.hpp"
#include "mvector.hpp"


/*! \brief Abstract base class for solving linear and nonlinear
 *  problems. Different implementation may exist.
 */
class Solver {

public:

    /*! \brief Virtual destructor.
     */
    virtual ~Solver() {}

    /*! \brief Solve problem \a p. Initial guess and solution are in
     *  vector \a X.
     */
    virtual void solve( const class Problem &p, Vector &X ) const = 0;

};


#endif













