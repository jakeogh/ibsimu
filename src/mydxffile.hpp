/*! \file mydxffile.hpp
 *  \brief Header file for mydxffile.hpp
 */

/* Copyright (c) 2010 Taneli Kalvas. All rights reserved.
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

#ifndef MY_DXF_FILE_HPP
#define MY_DXF_FILE_HPP 1



//#define MYDXF_DEBUG 1



#include <fstream>
#include "mydxfheader.hpp"
#include "mydxfentities.hpp"




/*! \brief DXF file class.
 *
 *  This class is a memory representation of a dxf file read from the
 *  disc.
 */
class MyDXFFile
{
    std::ifstream _fstr;
    bool _ascii;
    int  _linec;

    int          _group_code;

    int          _group_type;
    std::string  _group_string;
    double       _group_double;
    bool         _group_bool;
    int8_t       _group_int8;
    int16_t      _group_int16;
    int32_t      _group_int32;
    int64_t      _group_int64;

    class MyDXFHeader *_header;
    class MyDXFEntities *_entities;


public:
	
    MyDXFFile( const std::string &filename );
    ~MyDXFFile();

    int read_group( void );

    int group_get_code( void ) const;

    std::string group_get_string( void ) const;
    double group_get_double( void ) const;
    bool group_get_bool( void ) const;
    int8_t group_get_int8( void ) const;
    int16_t group_get_int16( void ) const;
    int32_t group_get_int32( void ) const;
    int64_t group_get_int64( void ) const;

    int linec( void ) const { return( _linec ); }

    class MyDXFEntities *get_entities( void ) { return( _entities ); };
};


#endif
