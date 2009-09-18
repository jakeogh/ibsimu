/*! \file bicgstab_solver.hpp
 *  \brief Header file for bicgstab_solver.hpp
 */

/* Copyright (c) 2005-2009 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * tvkalvas@cc.jyu.fi.
 * 
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
 */

#ifndef BICGSTAB_SOLVER_HPP
#define BICGSTAB_SOLVER_HPP 1


#include <iostream>
#include "solver.hpp"
#include "geometry.hpp"
#include "problem.hpp"
#include "matrix.hpp"
#include "mvector.hpp"


/*! \brief BiCGSTAB based solver implementation.
 *
 *  An implementation of virtual class Solver using ILU0
 *  preconditioner and bicgstab() solver function. Linear and
 *  nonlinear solvers.
 *
 */
class BiCGSTABSolver : public Solver {
    double   _eps;          /*!< \brief Accuracy request. */
    int      _imax;         /*!< \brief Maximum iteration count. */

    double   _newton_Reps;  /*!< \brief Accuracy request for Newton-Raphson residual. */
    double   _newton_dXeps; /*!< \brief Accuracy request for Newton-Raphson step. */
    double   _newton_imax;  /*!< \brief Maximum number of Newton-Raphson iterations. */

public:

    /*! \brief Constructor.
     */
    BiCGSTABSolver( double eps = 1.0e-6, int imax = 10000,
		    double newton_Reps = 1.0e-5, double newton_dXeps = 1.0e-6, int newton_imax = 10 );

    /*! \brief Destructor.
     */
    ~BiCGSTABSolver() {}

    /*! \brief Solve problem \a p defined in geometry \a g.  Initial
     *  guess and solution are in vector \a X.
     */
    virtual void solve( const Problem &p, Vector &X );

    /*! \brief Reset solver.
     *
     *  This is a signal from the problem that the problem has changed
     *  and internal caches (if they exist) in the solver should be
     *  resetted.
     */
    virtual void reset( void );

    /*! \brief Sets the accuracy request for BiCGSTAB solver.
     */
    void set_eps( double eps ) {
	_eps = eps;
    }

    /*! \brief Sets maximum iteration count for BiCGSTAB solver.
     */
    void set_imax( int imax ) {
	_imax = imax;
    }

    /*! \brief Sets maximum iteration count for Newton-Raphson steps.
     */
    void set_newton_imax( int newton_imax ) {
	_newton_imax = newton_imax;
    }

    /*! \brief Sets the accuracy request for Newton-Raphson residual.
     */
    void set_newton_residual_eps( int newton_Reps ) {
	_newton_Reps = newton_Reps;
    }

    /*! \brief Sets the accuracy request for Newton-Raphson step size.
     */
    void set_newton_step_eps( int newton_dXeps ) {
	_newton_dXeps = newton_dXeps;
    }

};


#endif















