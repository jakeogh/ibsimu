/*! \file epot_mgsolver.hpp
 *  \brief Multigrid solver for electric potential problem
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

#ifndef EPOT_MGSOLVER_HPP
#define EPOT_MGSOLVER_HPP 1


#include "epot_solver.hpp"


/*! \brief Multigrid solver for Electric potential problem.
 */
class EpotMGSolver : public EpotSolver {

    std::vector<MeshScalarField *>   _epotv;
    std::vector<Geometry *>          _geomv;
    std::vector<EpotSolver *>        _epotsolverv;
    std::vector<MeshScalarField *>   _rhsv;

    bool             _geom_prepared;  /*!< \brief Is geometry prepared? */
    uint32_t         _levels;         /*!< \brief Multigrid levels. */
    uint32_t         _npre;           /*!< \brief Pre cycle smoother rounds. */
    uint32_t         _npost;          /*!< \brief Post cycle smoother rounds. */
    uint32_t         _ncyc;           /*!< \brief Multigrid cycles. */
    double           _res;            /*!< \brief Residual error from top level. */
    double           _eps;            /*!< \brief Acceptable error for coarsest level. */
    
    
    double gs_loop_2d( void ) const;
    double gs_process_near_solid_2d( const uint8_t *nearsolid_ptr, 
				     uint32_t a, uint32_t dj ) const;
    double gs_process_pure_vacuum_2d( uint32_t a, uint32_t dj ) const;
    double gs_process_neumann_2d( uint32_t boundary, uint32_t a, uint32_t dj ) const;

    uint32_t number_of_dimensions( void ) const;
    void prepare_mg_geom( void );

    void preprocess( MeshScalarField &epot, const MeshScalarField &scharge );
    void postprocess( void );


    /*! \brief Reset solver/problem settings.
     */
    virtual void reset_problem( void );

    /*! \brief Solve problem with given mesh based space charge.
     */
    virtual void subsolve( MeshScalarField &epot, const MeshScalarField &scharge );

public:

    /*! \brief Constructor.
     */
    EpotMGSolver( Geometry &geom );

    /*! \brief Construct from file.
     */
    EpotMGSolver( Geometry &geom, std::istream &s );

    /*! \brief Destructor.
     */
    virtual ~EpotMGSolver();

    /*! \brief Sets the accuracy request.
     */
    void set_eps( double eps );

    /*! \brief Get estimate of residual error.
     */
    double get_residual( void ) const;

    /*! \brief Sets multigrid levels.
     */
    void set_levels( uint32_t levels );

    /*! \brief Sets number of multigrid cycles.
     */
    void set_ncyc( uint32_t ncyc );

    /*! \brief Sets number of pre cycle smoother rounds.
     */
    void set_npre( uint32_t npre );

    /*! \brief Sets number of post cycle smoother rounds.
     */
    void set_npost( uint32_t npost );

    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;

    /*! \brief Saves problem data to stream.
     */
    virtual void save( std::ostream &s ) const;
};

#endif
