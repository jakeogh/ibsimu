/*! \file file.hpp
 *  \brief Header file for file saving and loading functions
 */


#ifndef FILE_HPP
#define FILE_HPP 1


#include <iostream>
#include <zlib.h>
#include "error.hpp"


/* **************** *
 * Write            *
 * **************** */


/*! \brief Write int8_t \a value into stream \a s.
 */
inline void write_int8( std::ostream &s, int8_t value )
{
    s.write( (const char *)(&value), sizeof(int8_t) );
}


/*! \brief Write int16_t \a value into stream \a s.
 */
inline void write_int16( std::ostream &s, int16_t value )
{
    s.write( (const char *)(&value), sizeof(int16_t) );
}


/*! \brief Write int32_t \a value into stream \a s.
 */
inline void write_int32( std::ostream &s, int32_t value )
{
    s.write( (const char *)(&value), sizeof(int32_t) );
}


/*! \brief Write uint32_t \a value into stream \a s.
 */
inline void write_uint32( std::ostream &s, uint32_t value )
{
    s.write( (const char *)(&value), sizeof(uint32_t) );
}


/*! \brief Write double \a value into stream \a s.
 */
inline void write_double( std::ostream &s, double value )
{
    s.write( (const char *)(&value), sizeof(double) );
}


/*! \brief Write data block \a data of length \a len bytes into stream
 *  \a s in compressed form.
 */
inline void write_compressed_block( std::ostream &s, uint32_t len, const int8_t *data )
{
    uLongf  buflen = uLongf(len*1.1)+12 ;
    Bytef  *buf = new Bytef [buflen];
    compress( buf, &buflen, (Bytef *)data, (uLong)len );
    write_uint32( s, len );
    write_uint32( s, buflen );
    s.write( (const char *)buf, buflen );
    delete [] buf;
}



/* **************** *
 * Read             *
 * **************** */


inline int8_t read_int8( std::istream &s )
{
    int8_t value;
    s.read( (char *)(&value), sizeof(int8_t) );
    return( value );
}


inline int16_t read_int16( std::istream &s )
{
    int16_t value;
    s.read( (char *)(&value), sizeof(int16_t) );
    return( value );
}


inline int32_t read_int32( std::istream &s )
{
    int32_t value;
    s.read( (char *)(&value), sizeof(int32_t) );
    return( value );
}


inline uint32_t read_uint32( std::istream &s )
{
    uint32_t value;
    s.read( (char *)(&value), sizeof(uint32_t) );
    return( value );
}


inline double read_double( std::istream &s )
{
    double value;
    s.read( (char *)(&value), sizeof(double) );
    return( value );
}


inline uint32_t read_compressed_block( std::istream &s, uint32_t len, int8_t *dest )
{
    uint32_t datalen;
    uint32_t compressedlen;
    datalen = read_uint32( s );
    compressedlen = read_uint32( s );

    if( datalen > len )
	throw( Error( ERROR_LOCATION, "compressed data length longer than expected (" + 
		      to_string(datalen) + " > " + to_string(len) ) );

    Bytef  *buf = new Bytef [compressedlen];
    uLongf  tmp = datalen;
    s.read( (char *)buf, compressedlen );
    uncompress( (Bytef *)dest, &tmp, buf, compressedlen );
    delete [] buf;

    return( datalen );
}


#endif














