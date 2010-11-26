/*! \file mydxftables.hpp
 *  \brief Tables file for mydxftables.hpp
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

#ifndef MY_DXF_TABLES_HPP
#define MY_DXF_TABLES_HPP 1


#include <vector>
#include <string>
#include <stdint.h>
#include "mydxffile.hpp"


/*! \brief DXF table entry
 */
class MyDXFEntry
{

    std::string _handle;           // 5 (for others) or 105 (for DIMSTYLE)
    std::string _handle_to_owner;  // 330

protected:

    /*! \brief Constructor.
     */
    MyDXFEntry();

    /*! \brief Process group not belonging to the child entry.
     */
    void process_group( class MyDXFFile *dxf );
	
    /*! \brief Write common groups.
     */
    void write_common( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Debug print common groups.
     */
    void debug_print_common( std::ostream &os ) const;
    
public:

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFEntry() {};

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr ) = 0;

    /*! \brief Debug print.
     */
    virtual void debug_print( std::ostream &os ) const = 0;
    
    friend std::ostream &operator<<( std::ostream &os, const MyDXFEntry &e );
};



/*! \brief DXF table entry for block record table.
 */
class MyDXFEntry_BlockRecord : public MyDXFEntry
{

    std::string _name;             // 2
    int16_t     _units;            // 70
    int8_t      _explodability;    // 280
    int8_t      _scalability;      // 281
    std::string _handle_to_layout; // 340

public:

    /*! \brief Construct entry by reading from DXF file.
     */
    MyDXFEntry_BlockRecord( class MyDXFFile *dxf );

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFEntry_BlockRecord();

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Debug print.
     */
    virtual void debug_print( std::ostream &os ) const;
};



/*! \brief DXF table entry for layer table.
 */
class MyDXFEntry_Layer : public MyDXFEntry
{

    std::string _name;             // 2
    std::string _linetype;         // 6
    int16_t     _flags;            // 70
    int16_t     _color;            // 62 (negative if layer off)
    bool        _plotting;         // 290
    int8_t      _lineweight;       // 370

    std::string _handle_to_plot_style_name; // 390
    std::string _handle_to_material; // 347

public:

    /*! \brief Construct entry by reading from DXF file.
     */
    MyDXFEntry_Layer( class MyDXFFile *dxf );

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFEntry_Layer();

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Debug print.
     */
    virtual void debug_print( std::ostream &os ) const;
};



/*! \brief DXF table class.
 */
class MyDXFTable
{
    std::string               _name;            // 2
    std::string               _handle;          // 5
    std::string               _handle_to_owner; // 330

    std::vector<MyDXFEntry *> _entries;         // 70 (size)

public:

    MyDXFTable( const std::string &name, class MyDXFFile *dxf );
    ~MyDXFTable();

    /*! \brief Write dxf file to stream.
     */
    void write( class MyDXFFile *dxf, std::ofstream &ostr );

    void debug_print( std::ostream &os ) const;
};



/*! \brief DXF tables class.
 *
 *  Container for data of a DXF file tables.
 */
class MyDXFTables
{

    MyDXFTable  *_blockrecord;
    MyDXFTable  *_layer;

public:

    MyDXFTables( class MyDXFFile *dxf );
    ~MyDXFTables();

    /*! \brief Write dxf file to stream.
     */
    void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;
};





#endif



