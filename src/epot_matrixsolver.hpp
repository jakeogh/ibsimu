/*! \file epot_matrixsolver.hpp
 *  \brief Matrix solver for electric potential problem
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


#ifndef EPOT_MATRIXSOLVER_HPP
#define EPOT_MATRIXSOLVER_HPP 1


#include "epot_solver.hpp"
#include "crowmatrix.hpp"
#include "mvector.hpp"


/*! \brief Parent class for Matrix-based solvers for Electric potential problem.
 */
class EpotMatrixSolver : public EpotSolver {

protected:

    /*! \brief Class nodes to degrees of freedom mapping.
     *
     *  Uses running numbers starting from 0 to point to vectors with
     *  free variables (electric potential and matrices during solving
     *  for example). If a node is solid interior point or dirichlet
     *  boundary, a fixed negative number corresponding to the solid
     *  is used, i.e. from -1 to -6 for boundaries and starting from
     *  -7 for electrodes.
     */
    class Node2DoF {
        Int3D         _size;          /*!< \brief Size of mesh */
        int32_t      *_n2d;           /*!< \brief Nodes to degrees of freedom array. */
        
    public:
        
        Node2DoF();
        Node2DoF( Int3D size );
        ~Node2DoF();
	void clear( void );
        void resize( Int3D size );
	
        int32_t &operator()( int i ) 
            { return( _n2d[i] ); }
        int32_t &operator()( int i, int j ) 
            { return( _n2d[i+j*_size[0]] ); }
        int32_t &operator()( int i, int j, int k ) 
            { return( _n2d[i+j*_size[0]+k*_size[0]*_size[1]] ); }
        
        const int32_t &operator()( int i ) const
            { return( _n2d[i] ); }
        const int32_t &operator()( int i, int j )  const 
            { return( _n2d[i+j*_size[0]] ); }
        const int32_t &operator()( int i, int j, int k ) const 
            { return( _n2d[i+j*_size[0]+k*_size[0]*_size[1]] ); }
        
        /*! \brief Print debugging information to os.
         */
        void debug_print( std::ostream &os ) const;
    };

    uint32_t               _dof;           /*!< \brief Degrees of freedom. */
    Node2DoF               _n2d;           /*!< \brief Nodes to degrees of freedom map. */
    CRowMatrix            *_fd_mat;        /*!< \brief Finite Difference matrix. */
    Vector                *_fd_vec;        /*!< \brief Finite Difference vector. */

    MeshScalarField       *_epot;
    const MeshScalarField *_scharge;
    

    /*! \brief Constructor.
     */
    EpotMatrixSolver( Geometry &geom );

    /*! \brief Construct from file.
     */
    EpotMatrixSolver( Geometry &geom, std::istream &s );

    /*! \brief Return const pointers to the matrix \a A and vector \a
     *  B of the linear problem.
     *
     *  This function should only be called by solver, while method
     *  subsolve() is running and only in case of linear problem.
     */
    void get_vecmat( const Matrix **A, const Vector **B ) const;

    /*! \brief Return const pointers to jacobian matrix and residual
     *  vector of the problem to \a J and \a R at \a X.
     *
     *  This function should only be called by solver, while method
     *  subsolve() is running and only in case of nonlinear problem.
     */
    void get_resjac( const Matrix **J, const Vector **R, const Vector &X ) const;

    /*! \brief Return true if problem is linear.
     */
    bool linear( void ) const;

    /*! \brief Preprocess.
     *
     *  Modify solid mesh suitable for solver, build n2d map and
     *  linear matrix. Make right-hand-side and solution vector from
     *  electric potential.
     */
    void preprocess( MeshScalarField &epot, const MeshScalarField &scharge, Vector &X );

    /*! \brief Postprocess.
     *
     *  Return solid mesh back to original state and make electric
     *  potential field from solution vector.
     */
    void postprocess( MeshScalarField &epot, const Vector &X );

    /*! \brief Reset matrix representation.
     */
    void reset_matrix( void );

private:

    void set_link( int32_t a, int32_t b, double val );
    void add_vacuum_node( uint32_t i, uint32_t j, uint32_t k );
    void add_near_solid_node_1d( uint32_t i );
    void add_near_solid_node_2d( uint32_t i, uint32_t j );
    void add_near_solid_node_cyl( uint32_t i, uint32_t j );
    void add_near_solid_node_3d( uint32_t i, uint32_t j, uint32_t k );
    void add_near_solid_node( uint32_t i, uint32_t j, uint32_t k );
    void add_neumann_node( uint32_t i, uint32_t j, uint32_t k, uint32_t boundary );

public:

    /*! \brief Destructor.
     */
    virtual ~EpotMatrixSolver();





    /*! \brief Print debugging information to os.
     */
    virtual void debug_print( std::ostream &os ) const;

    /*! \brief Saves problem data to stream.
     */
    virtual void save( std::ostream &s ) const;
};


#endif
