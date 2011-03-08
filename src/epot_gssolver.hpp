/*! \file epot_gssolver.hpp
 *  \brief Gauss-Seidel solver for electric potential problem
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

#ifndef EPOT_GSSOLVER_HPP
#define EPOT_GSSOLVER_HPP 1


#include "epot_solver.hpp"


class EpotGSSolver : public EpotSolver {

    MeshScalarField *_epot;
    MeshScalarField *_rhs;

    uint32_t         _iter;           /*!< \brief Number of iteration rounds done. */
    uint32_t         _imax;           /*!< \brief Maximum number of iteration rounds. */
    double           _eps;            /*!< \brief Accuracy request. */
    double           _res;            /*!< \brief Residual error. */
    double           _w;              /*!< \brief Relaxation coefficient. */

    std::vector<uint32_t> _nsind;     /*!< \brief Stored near solid indexes for Neumann conversion nodes. */

    double gs_loop_3d( void ) const;
    double gs_process_near_solid_3d( const uint8_t *nearsolid_ptr, 
				     uint32_t a, uint32_t dj, uint32_t dk ) const;
    double gs_process_pure_vacuum_3d( uint32_t a, uint32_t dj, uint32_t dk ) const;
    double gs_process_neumann_3d( uint32_t boundary, uint32_t a,
				  uint32_t dj, uint32_t dk ) const;


    double gs_loop_cyl( void ) const;
    double gs_process_near_solid_cyl( const uint8_t *nearsolid_ptr, 
				     uint32_t i, uint32_t j ) const;
    double gs_process_pure_vacuum_cyl( uint32_t i, uint32_t j ) const;
    double gs_process_neumann_cyl( uint32_t boundary, uint32_t a, uint32_t dj ) const;


    double gs_loop_2d( void ) const;
    double gs_process_near_solid_2d( const uint8_t *nearsolid_ptr, 
				     uint32_t a, uint32_t dj ) const;
    double gs_process_pure_vacuum_2d( uint32_t a, uint32_t dj ) const;
    double gs_process_neumann_2d( uint32_t boundary, uint32_t a, uint32_t dj ) const;


    double gs_loop_1d( void ) const;
    double gs_process_near_solid_1d( const uint8_t *nearsolid_ptr, 
				     uint32_t i ) const;
    double gs_process_pure_vacuum_1d( uint32_t i ) const;
    double gs_process_neumann_1d( uint32_t boundary, uint32_t i ) const;


    void preprocess( const MeshScalarField &scharge );
    void postprocess( void );
    virtual void subsolve( MeshScalarField &epot, const MeshScalarField &scharge );

public:

    /*! \brief Constructor.
     */
    EpotGSSolver( Geometry &geom );

    /*! \brief Construct from file.
     */
    EpotGSSolver( Geometry &geom, std::istream &s );

    /*! \brief Destructor.
     */
    virtual ~EpotGSSolver() {}

    /*! \brief Sets the accuracy request.
     */
    void set_eps( double eps );

    /*! \brief Get estimate of residual error.
     */
    double get_residual( void ) const;

    /*! \brief Sets maximum iteration count.
     */
    void set_imax( uint32_t imax );

    /*! \brief Sets relaxation parameter.
     */
    void set_w( double w );

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;

    /*! \brief Saves problem data to stream.
     */
    virtual void save( std::ostream &s ) const;
};

#endif
