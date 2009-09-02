/*! \file matrix.hpp
 *  \brief Header file for matrix class.
 */


#ifndef MATRIX_HPP
#define MATRIX_HPP 1


#include <iostream>
#include "mvector.hpp"


/*! \brief Container object for matrix-vector multiplication operation.
 *
 *  This container object is used to store a matrix-vector
 *  multiplication operation. For more information about the use of
 *  MCRowMatrixMulVec, see MCRowMatrix.
 */
struct MatrixMulVec {
    const class Matrix  *_mat; //!< Pointer to matrix.
    const class Vector  *_vec; //!< Pointer to vector.

    /*! \brief Constructor for MMatrixMulVec with matrix \a mat and
     *  vector \a vec.
     */
    MatrixMulVec( const Matrix &mat, const class Vector &vec ) : 
	_mat(&mat), _vec(&vec) {}

    friend class Vector;
};


/*! \brief Base matrix class.
 *
 *  The matrix class is an abstract class designed to be used as a
 *  base class for different matrix implementations.
 */
class Matrix {
    virtual double get_check( uint32_t i, uint32_t j ) const = 0;
    virtual double &set_check( uint32_t i, uint32_t j ) = 0;
    virtual double get_no_check( uint32_t i, uint32_t j ) const = 0;
    virtual double &set_no_check( uint32_t i, uint32_t j ) = 0;

public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    /*! \brief Virtual destructor.
     */
    virtual ~Matrix() {}

/* ************************************** *
 * Access and information                 *
 * ************************************** */

    /*! \brief Returns the number of columns of the matrix.
     */
    virtual uint32_t columns( void ) const = 0;

    /*! \brief Returns the number of rows of the matrix.
     */
    virtual uint32_t rows( void ) const = 0;

    /*! \brief Returns the number of rows \a n and the number of
     *  columns \a m of the matrix.
     */
    virtual void size( uint32_t &n, uint32_t &m ) const = 0;

/* ************************************** *
 * User level control                     *
 * ************************************** */

    /*! \brief Resizes the matrix to \a nn x \a mm.
     */
    virtual void resize( uint32_t n, uint32_t m ) = 0;

    //virtual void merge( Matrix &mat ) = 0;

    /*! \brief Clears the matrix (sets all element to zero).
     */
    virtual void clear( void ) = 0;

/* ************************************** *
 * User level matrix element access       *
 * ************************************** */

    /*! \brief Function to get a matrix element value at (i,j).
     */
    inline double get( uint32_t i, uint32_t j ) const;

    /*! \brief Function to get a reference to matrix element value at
     *  (i,j).
     */
    inline double &set( uint32_t i, uint32_t j );

/* ************************************** *
 * Matrix-Vector operations               *
 * ************************************** */

    /*! \brief Operator for matrix-vector multiplication.
     */
    MatrixMulVec operator*( const class Vector &vec ) const;

    virtual void multiply_by_vector( Vector &res, const Vector &rhs ) const = 0;
    virtual void lower_unit_solve( Vector &y, const Vector &b ) const = 0;
    virtual void upper_diag_solve( Vector &x, const Vector &y ) const = 0;

    friend class Vector;
};


inline double Matrix::get( uint32_t i, uint32_t j ) const
{
#ifdef SPM_RANGE_CHECK
    return( get_check( i, j ) );
#else
    return( get_no_check( i, j ) );
#endif
}    


inline double &Matrix::set( uint32_t i, uint32_t j )
{
#ifdef SPM_RANGE_CHECK
    return( set_check( i, j ) );
#else
    return( set_no_check( i, j ) );
#endif
}    


std::ostream &operator<<( std::ostream &os, const Matrix &mat );



#endif














