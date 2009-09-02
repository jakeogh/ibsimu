/*! \file interpolation.hpp
 *  \brief Header file for 2d interpolations.
 */


#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP 1


#include <vector>


/*! \brief Base class for 2d interpolation.
 *
 *  Provides an interpolation of a function defined at mesh points \a
 *  (x,y), where 0 <= (x,y) <= 1. 
 */
class Interpolation2D
{
protected:

    size_t                    _n;    /*!< \brief Size of first coordinate of mesh. */
    size_t                    _m;    /*!< \brief Size of second coordinate of mesh. */
    std::vector<double>       _f;    /*!< \brief Function data mesh. */

    /*! \brief Constructor.
     */
    Interpolation2D( size_t n, size_t m, const std::vector<double> &f );

    const double &__f( int i, int j ) const;
    double &__f( int i, int j );

public:

    /*! \brief Virtual destructor.
     */
    virtual ~Interpolation2D() {}

    /*! \brief Operator for getting interpolation at \a (x,y).
     *
     *  Returns an interpolated value of the function at \a (x,y),
     *  where 0 <= (x,y) <= 1. Returns NaN outside this area.
     */
    virtual double operator()( double x, double y ) const = 0;
};





/*! \brief Closest point 2d interpolation.
 *
 *  Not really an interpolation. Just returns the closest point of
 *  original data.
 */
class ClosestInterpolation2D : public Interpolation2D
{

public:

    /*! \brief Constructor.
     */
    ClosestInterpolation2D( size_t n, size_t m, const std::vector<double> &f );

    /*! \brief Destructor.
     */
    virtual ~ClosestInterpolation2D() {}

    /*! \brief Operator for getting interpolation at \a (x,y).
     *
     *  Returns an interpolated value of the function at \a (x,y),
     *  where 0 <= (x,y) <= 1. Returns NaN outside this area.
     */
    virtual double operator()( double x, double y ) const;
};


/*! \brief BiLinear 2d interpolation.
 */
class BiLinearInterpolation2D : public Interpolation2D
{

public:

    /*! \brief Constructor.
     */
    BiLinearInterpolation2D( size_t n, size_t m, const std::vector<double> &f );

    /*! \brief Destructor.
     */
    virtual ~BiLinearInterpolation2D() {}

    /*! \brief Operator for getting interpolation at \a (x,y).
     *
     *  Returns an interpolated value of the function at \a (x,y),
     *  where 0 <= (x,y) <= 1. Returns NaN outside this area.
     */
    virtual double operator()( double x, double y ) const;
};





/*! \brief BiCubic 2d interpolation.
 *
 *  Calculates the derivatives of the function at mesh points using
 *  central finite differences. Zero derivatives are assumed at boundaries.
 */
class BiCubicInterpolation2D : public Interpolation2D
{

    std::vector<double>       _fx;   /*!< \brief X-derivative of function, data mesh. */
    std::vector<double>       _fy;   /*!< \brief Y-derivative of function, data mesh. */
    std::vector<double>       _fxy;  /*!< \brief XY-derivative of function, data mesh. */
    std::vector<double>       _c;    /*!< \brief Precalculated coefficients for interpolation. 
				     *
				     * 16 numbers per mesh square. Totalling 16*(n-1)*(m-1).
				     */

    const double &__fx( int i, int j ) const;
    const double &__fy( int i, int j ) const;
    const double &__fxy( int i, int j ) const;

    double &__fx( int i, int j );
    double &__fy( int i, int j );
    double &__fxy( int i, int j );

    static void calc_coefs( double *c, double *x );
    static const double wt[16][16];

public:

    /*! \brief Constructor.
     */
    BiCubicInterpolation2D( size_t n, size_t m, const std::vector<double> &f );

    /*! \brief Destructor.
     */
    virtual ~BiCubicInterpolation2D() {}

    /*! \brief Operator for getting interpolation at \a (x,y).
     *
     *  Returns an interpolated value of the function at \a (x,y),
     *  where 0 <= (x,y) <= 1. Returns NaN outside this area.
     */
    virtual double operator()( double x, double y ) const;
};


#endif














