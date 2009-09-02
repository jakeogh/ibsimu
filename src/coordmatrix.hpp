/*! \file coordmatrix.hpp
 *  \brief Header file for the coordinate sparse matrix class.
 */


#ifndef COORDMATRIX_HPP
#define COORDMATRIX_HPP 1


#include <cstdlib>
#include <iostream>
#include "matrix.hpp"
#include "error.hpp"


/*! \brief Coordinate sparse matrix class.
 *
 *  The matrix is stored in the standard coordinate sparse matrix
 *  storage mode.  In coordinate storage method all non-zero
 *  matrix elements are stored in array as triplets 
 *  (\a row, \a column, \a val). The format itself does not require
 *  a certain ordering of elements, but some implementation might 
 *  need/be faster on some ordering. Our example matrix \code
 *      | 1  2  0  0  3|
 *      | 4  5  6  0  0|
 *  A = | 0  7  8  0  9|
 *      | 0  0  0 10  0|
 *      |11  0  0  0 12|
 *  \endcode is represented in compressed column sparse matrix class as: \code
 *  val[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}
 *  row[] = {0, 0, 0, 1, 1, 1, 2, 2, 2,  3,  4,  4}
 *  col[] = {0, 1, 4, 0, 1, 2, 1, 2, 4,  3,  0,  4}
 *  \endcode
 */
class CoordMatrix : public Matrix {
    uint32_t  _n;      //!< Number of rows.
    uint32_t  _m;      //!< Number of columnss.
    uint32_t  _nz;     //!< Number of nonzero elements.
    uint32_t  _asize;  //!< Allocation size of \a _col and \a _val.
    uint32_t *_row;    //!< Row indices(i), \a _nz elements in use, \a _asize elements allocated.
    uint32_t *_col;    //!< Column indices(j), \a _nz elements in use, \a _asize elements allocated.
    double   *_val;    //!< Element values, \a _nz elements in use, \a _asize elements allocated.

    void allocate( void );
    void reallocate( void );

    double get_check( uint32_t i, uint32_t j ) const;
    double &set_check( uint32_t i, uint32_t j );
    double get_no_check( uint32_t i, uint32_t j ) const;
    double &set_no_check( uint32_t i, uint32_t j );

    void clear_check( uint32_t i, uint32_t j );
    void clear_no_check( uint32_t i, uint32_t j );

public:

/* ************************************** *
 * Constructors and destructor            *
 * ************************************** */

    CoordMatrix() : _n(0), _m(0), _nz(0), _asize(0), _row(NULL), _col(NULL), _val(NULL) { }
    CoordMatrix( uint32_t n, uint32_t m );
    CoordMatrix( uint32_t n, uint32_t m, uint32_t nz, 
		 const uint32_t *row, const uint32_t *col, const uint32_t *val );
    CoordMatrix( const CoordMatrix &mat );
    CoordMatrix( const class CRowMatrix &mat );
    CoordMatrix( const class CColMatrix &mat );
    ~CoordMatrix();

/* ************************************** *
 * Access and information                 *
 * ************************************** */

    /*! \brief Returns the number of columns in the matrix.
     */
    uint32_t columns( void ) const { return( _m ); }

    /*! \brief Returns the number of rows in the matrix.
     */
    uint32_t rows( void ) const { return( _n ); }

    /*! \brief Returns the number of columns and number of columns in \a nn and \a mm.
     */
    void size( uint32_t &n, uint32_t &m ) const { n = _n; m = _m; }

    /*! \brief Returns the number of non-zero elements in the matrix.
     */
    uint32_t nz_elements( void ) const { return( _nz ); }

    /*! \brief Returns the number of elements allocated for matrix.
     */
    uint32_t capacity( void ) const { return( _asize ); }

/* ************************************** *
 * User level control                     *
 * ************************************** */

    /*! \brief Resizes the matrix to size \a n x \a m.
     *
     *  All existing non-zero elements are cleared.
     */
    void resize( uint32_t n, uint32_t m );

    /*! \brief Merges matrix \a mat into the matrix leaving \a mat empty.
     *
     *  Copies contents of matrix \a mat into the matrix and sets
     *  contents of matrix \a mat to \a n = 0 and \a m = 0.
     *  \param mat Matrix to copy from.
     */
    void merge( CoordMatrix &mat );

    /*! \brief Clear non-zero matrix elements, set all elements to zero.
     */
    void clear( void );

    /*! \brief Clear matrix element (i,j).
     *
     *  Removes element (i,j) from the list of non-zero matrix elements.
     */
    void clear( uint32_t i, uint32_t j );

    /*! \brief Reserve memory for \a size matrix elements.
     */
    void reserve( uint32_t size );

    /*! \brief Order (sort) matrix data in ascending (row,column)
     *  index order.
     */
    void order_ascending_row_column( void );

    /*! \brief Order (sort) matrix data in ascending (column,row)
     *  index order.
     */
    void order_ascending_column_row( void );

    /*! \brief Prints the values of all internal data to std::cout.
     */
    void debug_print( void ) const;

/* ************************************** *
 * User level matrix element access       *
 * ************************************** */

    /*! \brief Function to get a matrix element value at (i,j).
     *
     *  Range checking is done for \a i and \a j if \c SPM_RANGE_CHECK
     *  is defined. Throws ErrorRange exception on range checking errors.
     */
    double get( uint32_t i, uint32_t j ) const;

    /*! \brief Function to get a reference to matrix element value at (i,j).
     *
     *  This function can be used to set or modify matrix element
     *  value. See following examples: \code
     *  A.set(0,0) = 1.2
     *  A.set(0,1) *= 2
     *  A.set(0,1) += 0.1 
     *  \endcode Note that a
     *  reference is actually a pointer to the memory location of the
     *  element and therefore unexpected things can happen if matrix
     *  is modified while using set, for example \code 
     *  A.set(0,0) = A.set(0,1) = A.set(0,2) = 5.0 
     *  \endcode does not do what you
     *  would expect it to do.
     *
     *  Range checking is done for \a i and \a j if \c SPM_RANGE_CHECK
     *  is defined. Throws ErrorRange exception on range checking errors.
     */
    double &set( uint32_t i, uint32_t j );

    /*! \brief Set element with no checks.
     *
     * Sets element (\a i, \a j) to value \a vval. No checks are
     * performed. User must be certain that duplicate entries are not
     * made.
     */
    void set_no_duplicate_check( uint32_t i, uint32_t j, double vval );

/* ************************************** *
 * Low level access                       *
 * ************************************** */

    /*! \brief Returns a reference to the to the internal row data
     *  of the matrix.
     */
    uint32_t &row( uint32_t i ) { return( _row[i] ); }

    /*! \brief Returns a reference to the to the internal column data
     *  \a ptr of the matrix.
     */
    uint32_t &col( uint32_t i ) { return( _col[i] ); }

    /*! \brief Returns a reference to the to the internal value data
     *  of the matrix.
     */
    double &val( uint32_t i ) { return( _val[i] ); }

    /*! \brief Returns a const reference to the to the internal row
     *  data of the matrix.
     */
    const uint32_t &row( uint32_t i ) const { return( _row[i] ); }

    /*! \brief Returns a const reference to the to the internal column
     *  data \a ptr of the matrix.
     */
    const uint32_t &col( uint32_t i ) const { return( _col[i] ); }

    /*! \brief Returns a const reference to the to the internal value
     *  data of the matrix.
     */
    const double &val( uint32_t i ) const { return( _val[i] ); }

    /*! \brief Set number of non-zero elements in the matrix.
     *
     *  This function is to be used with low level access
     *  functions. Internal arrays are resized if \a nz is larger than
     *  the allocated size.
     */
    void set_nz( uint32_t nz );

/* ************************************** *
 * Assignent operators                    *
 * ************************************** */

    CoordMatrix &operator=( const CoordMatrix &mat );
    CoordMatrix &operator=( const CColMatrix &mat );
    CoordMatrix &operator=( const CRowMatrix &mat );

/* ************************************** *
 * Matrix-Vector operations               *
 * ************************************** */

    /*  \brief Calculates \a x = \a A*b.
     */
    void multiply_by_vector( Vector &res, const Vector &rhs ) const;
    void lower_unit_solve( Vector &y, const Vector &b ) const;
    void upper_diag_solve( Vector &x, const Vector &y ) const;


    friend class CRowMatrix;
    friend class CColMatrix;
};


inline double CoordMatrix::get( uint32_t i, uint32_t j ) const
{
#ifdef SPM_RANGE_CHECK
    return( get_check( i, j ) );
#else
    return( get_no_check( i, j ) );
#endif
}    


inline double &CoordMatrix::set( uint32_t i, uint32_t j )
{
#ifdef SPM_RANGE_CHECK
    return( set_check( i, j ) );
#else
    return( set_no_check( i, j ) );
#endif
}    


inline void CoordMatrix::clear( uint32_t i, uint32_t j )
{
#ifdef SPM_RANGE_CHECK
    clear_check( i, j );
#else
    clear_no_check( i, j );
#endif
}


#endif
















