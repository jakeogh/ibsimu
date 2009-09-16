/*! \file vec3d.hpp
 *  \brief Header file defining classes Vec3D and Int3D.
 */


#ifndef VEC3D_HPP
#define VEC3D_HPP 1


#include <math.h>
#include <stdint.h>
#include <iostream>
#include <iostream>
#include <iomanip>
#include "file.hpp"


/*! \brief 3D %Vector class.
 */
class Vec3D {
    double p[3];

public:

    Vec3D() { p[0] = 0.0; p[1] = 0.0; p[2] = 0.0; }
    Vec3D( double x ) { p[0] = x; p[1] = 0.0; p[2] = 0.0; }
    Vec3D( double x, double y ) { p[0] = x; p[1] = y; p[2] = 0.0; }
    Vec3D( double x, double y, double z ) { p[0] = x; p[1] = y; p[2] = z; }
    Vec3D( std::istream &s ) {
	p[0] = read_double( s );
	p[1] = read_double( s );
	p[2] = read_double( s );
    }
    ~Vec3D() {}

    double &operator[]( int i ) { return( p[i] ); }
    const double &operator[]( int i ) const { return( p[i] ); }
    double &operator()( int i ) { return( p[i] ); }
    const double &operator()( int i ) const { return( p[i] ); }

    /*! \brief %Vector addition
     */
    Vec3D operator+( const Vec3D &vec ) const { 
	return( Vec3D( p[0] + vec[0], 
		       p[1] + vec[1],
		       p[2] + vec[2] ) );
    }

    /*! \brief %Vector difference
     */
    Vec3D operator-( const Vec3D &vec ) const {
	return( Vec3D( p[0] - vec[0],
		       p[1] - vec[1],
		       p[2] - vec[2] ) );
    } 

    /*! \brief %Vector accumulation
     */
    Vec3D &operator+=( const Vec3D &vec ) { 
	p[0] += vec[0];
	p[1] += vec[1];
	p[2] += vec[2];
	return( *this );
    }

    /*! \brief Dot product
     */
    double operator*( const Vec3D &vec ) const { 
	return( p[0] * vec[0] +
		p[1] * vec[1] +
		p[2] * vec[2] );
    }

    Vec3D operator*( double x ) { 
	return( Vec3D( x*p[0], x*p[1], x*p[2] ) );
    }

    Vec3D &operator*=( double x ) { 
	p[0] *= x;
	p[1] *= x;
	p[2] *= x;
	return( *this );
    }

    Vec3D &operator/=( double x ) { 
	double div = 1.0/x;
	p[0] *= div;
	p[1] *= div;
	p[2] *= div;
	return( *this );
    }

    bool operator!=( const Vec3D &x ) { 
	if( p[0] != x.p[0] || p[1] != x.p[1] || p[2] != x.p[2] )
	    return( true );
	return( false ); 
    }

    /*! \brief Normalize vector
     */
    void normalize() {
	double inv_norm = 1.0/sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
	p[0] *= inv_norm;
	p[1] *= inv_norm;
	p[2] *= inv_norm;
    }

    /*! \brief Returns 2-norm of vector
     *
     *  \f$ ||x||_2 = \sqrt{ \Sigma_{i=1}^n x_i^2 } \f$
     */
    double norm2() const {
	return( sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] ) );
    }

    /*! \brief Returns square of 2-norm of vector
     *
     *  \f$ (||x||_2)^2 = \Sigma_{i=1}^n x_i^2 \f$
     */
    double ssqr() const {
	return( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
    }

    void save( std::ostream &s ) const { 
	write_double( s, p[0] );
	write_double( s, p[1] );
	write_double( s, p[2] ); 
    }

    /*! \brief Cross product
     */
    friend Vec3D cross( const Vec3D &vec1, const Vec3D &vec2 );

    friend double norm2( const Vec3D &vec );

    friend Vec3D operator*( double x, const Vec3D &vec );
    friend Vec3D operator*( double x, const class Int3D &i );
    friend std::ostream &operator<<( std::ostream &os, const Vec3D &vec );
};


inline double norm2( const Vec3D &vec ) {
    return( vec.norm2() );
}

inline Vec3D cross( const Vec3D &vec1, const Vec3D &vec2 ) { 
    return( Vec3D( vec1[1] * vec2[2] - vec1[2] * vec2[1], 
		   vec1[2] * vec2[0] - vec1[0] * vec2[2],
		   vec1[0] * vec2[1] - vec1[1] * vec2[0] ) );
}


inline Vec3D operator*( double x, const Vec3D &vec )
{
    return( Vec3D( x*vec.p[0], x*vec.p[1], x*vec.p[2] ) );
}


inline std::ostream &operator<<( std::ostream &os, const Vec3D &vec ) 
{
    os << std::setw(12) << to_string(vec[0]).substr(0,12) << " ";
    os << std::setw(12) << to_string(vec[1]).substr(0,12) << " ";
    os << std::setw(12) << to_string(vec[2]).substr(0,12);
    return( os );
}


/*! \brief 3D Integer vector class.
 */
class Int3D {
    int32_t l[3];

public:

    Int3D() { l[0] = 0; l[1] = 0; l[2] = 0; }
    Int3D( int32_t i ) { l[0] = i; l[1] = 0; l[2] = 0; }
    Int3D( int32_t i, int32_t j ) { l[0] = i; l[1] = j; l[2] = 0; }
    Int3D( int32_t i, int32_t j, int32_t k ) { l[0] = i; l[1] = j; l[2] = k; }
    Int3D( std::istream &s ) {
	l[0] = read_int32( s );
	l[1] = read_int32( s );
	l[2] = read_int32( s );
    }
    ~Int3D() {}

    int32_t &operator[]( int i ) { return( l[i] ); }
    const int32_t &operator[]( int i ) const { return( l[i] ); }
    int32_t &operator()( int i ) { return( l[i] ); }
    const int32_t &operator()( int i ) const { return( l[i] ); }

    Int3D operator-( const Int3D &i ) {
	return( Int3D( l[0] - i[0],
		       l[1] - i[1],
		       l[2] - i[2] ) );
    } 

    Vec3D operator*( double x ) { 
	return( Vec3D( x*l[0], x*l[1], x*l[2] ) );
    }

    bool operator!=( const Int3D &i ) { 
	if( l[0] != i.l[0] || l[1] != i.l[1] || l[2] != i.l[2] )
	    return( true );
	return( false ); 
    }

    void save( std::ostream &s ) const { 
	write_int32( s, l[0] );
	write_int32( s, l[1] );
	write_int32( s, l[2] ); 
    }

    friend Vec3D operator*( double x, const Int3D &i );
    friend std::ostream &operator<<( std::ostream &os, const Vec3D &vec );
};


inline Vec3D operator*( double x, const Int3D &i )
{
    Vec3D res;
    res[0] = x*i.l[0];
    res[1] = x*i.l[1];
    res[2] = x*i.l[2];
    return( res );
}


inline std::ostream &operator<<( std::ostream &os, const Int3D &vec ) 
{
    os << std::setw(12) << to_string(vec[0]).substr(0,12) << " ";
    os << std::setw(12) << to_string(vec[1]).substr(0,12) << " ";
    os << std::setw(12) << to_string(vec[2]).substr(0,12);
    return( os );
}


#endif













