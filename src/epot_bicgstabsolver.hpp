/*! \file epot_bicgstabsolver.hpp
 *  \brief BiCGSTAB matrix solver for electric potential problem
 */

/* Copyright (c) 2011 Taneli Kalvas. All rights reserved.
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
 * taneli.kalvas@jyu.fi.
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


#ifndef EPOT_BICGSTABSOLVER_HPP
#define EPOT_BICGSTABSOLVER_HPP 1


#include "epot_matrixsolver.hpp"
#include "ccolmatrix.hpp"


/*! \brief BiCGSTAB matrix solver for Electric potential problem.
 */
class EpotBiCGSTABSolver : public EpotMatrixSolver {

    double   _eps;          /*!< \brief Accuracy request. */
    uint32_t _imax;         /*!< \brief Maximum iteration count. */

    uint32_t _iter;         /*!< \brief Number of iteration rounds done. */
    double   _res;          /*!< \brief Residual error. */

    double   _newton_Reps;  /*!< \brief Accuracy request for Newton-Raphson residual. */
    double   _newton_dXeps; /*!< \brief Accuracy request for Newton-Raphson step. */
    uint32_t _newton_imax;  /*!< \brief Maximum number of Newton-Raphson iterations. */

    /*! \brief Reset solver/problem settings.
     */
    virtual void reset_problem( void );

    /*! \brief Solve problem with given mesh based space charge.
     */
    virtual void subsolve( MeshScalarField &epot, const MeshScalarField &scharge );

public:

    /*! \brief Constructor.
     */
    EpotBiCGSTABSolver( Geometry &geom,
			double eps = 1.0e-4, 
			uint32_t imax = 10000,
			double newton_Reps = 1.0e-4, 
			double newton_dXeps = 1.0e-6, 
			uint32_t newton_imax = 10 );

    /*! \brief Construct from file.
     */
    EpotBiCGSTABSolver( Geometry &geom, std::istream &s );

    /*! \brief Destructor.
     */
    virtual ~EpotBiCGSTABSolver();

    /*! \brief Sets the accuracy request for BiCGSTAB solver.
     */
    void set_eps( double eps );

    /*! \brief Sets maximum iteration count for BiCGSTAB solver.
     */
    void set_imax( uint32_t imax );

    /*! \brief Sets maximum iteration count for Newton-Raphson steps.
     */
    void set_newton_imax( uint32_t newton_imax );

    /*! \brief Sets the accuracy request for Newton-Raphson residual.
     */
    void set_newton_residual_eps( double newton_Reps );

    /*! \brief Sets the accuracy request for Newton-Raphson step size.
     */
    void set_newton_step_eps( double newton_dXeps );

    /*! \brief Get estimate of residual error.
     */
    double get_residual( void ) const;

    /*! \brief Get number of iteration rounds done with last solve().
     */
    uint32_t get_iter( void ) const;

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;

    /*! \brief Saves problem data to stream.
     */
    virtual void save( std::ostream &s ) const;
};


#endif
