/*! \file mydxfentities.hpp
 *  \brief Header file for mydxfentities.hpp
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

#ifndef MY_DXF_ENTITIES_HPP
#define MY_DXF_ENTITIES_HPP 1


#include <stdint.h>
#include <vector>
#include "mydxffile.hpp"
#include "mydxfvec.hpp"


/*! \brief DXF entity base class.
 *
 *  A general base class for all DXF entities.
 */
class MyDXFEntity
{

protected: 

    std::string _handle;
    std::string _layer;

    MyDXFEntity();

    void process_group( class MyDXFFile *dxf );
    void constructor_debug_print( void ) const;
    void debug_print_base( std::ostream &os ) const;
    virtual void debug_print( std::ostream &os ) const = 0;
    
public:

    virtual ~MyDXFEntity() {}

    virtual MyDXFEntity *copy( void ) const = 0;

    void set_layer( const std::string &layer ) { _layer = layer; }
    std::string get_layer( void ) const { return( _layer ); }

    virtual void scale( double s ) = 0;

    friend std::ostream &operator<<( std::ostream &os, const MyDXFEntity &ent );
};


/*! \brief DXF path entity base class.
 *
 *  A base class for two dimensional DXF entities, which can
 *  be part of a path.
 */
class MyDXFPathEntity : public MyDXFEntity
{

protected: 

    MyDXFPathEntity() {}

public:

    virtual ~MyDXFPathEntity() {}

    virtual MyDXFVec start( void ) const = 0;
    virtual MyDXFVec end( void ) const = 0;

    virtual void set_start( const MyDXFVec &s ) = 0;
    virtual void set_end( const MyDXFVec &e ) = 0;

    /*! \brief Check for ray crossing.
     *
     *  Check if ray going from point (x,y) downwards (negative y
     *  direction) crosses the entity. Return 1 if crosses odd number
     *  of times and 0 if even number of times. Return 2 in case of
     *  exact crossing at boundaries. This function is used as a
     *  subroutine to inside_loop().
     */
    virtual int ray_cross( double x, double y ) const = 0;
};


/*! \brief DXF line entity class.
 *
 *
 */
class MyDXFLine : public MyDXFPathEntity
{

    MyDXFVec _p1;
    MyDXFVec _p2;

    virtual void debug_print( std::ostream &os ) const;

public:

    MyDXFLine( class MyDXFFile *dxf );
    virtual ~MyDXFLine() {}

    virtual MyDXFLine *copy( void ) const { return( new MyDXFLine( *this ) ); }

    virtual MyDXFVec start( void ) const { return( _p1 ); }
    virtual MyDXFVec end( void ) const { return( _p2 ); }

    virtual void set_start( const MyDXFVec &s );
    virtual void set_end( const MyDXFVec &e );

    /*! \brief Check for ray crossing.
     *
     *  Check if ray going from point (x,y) downwards (negative y
     *  direction) crosses the entity. Return 1 if crosses odd number
     *  of times and 0 if even number of times. Return 2 in case of
     *  exact crossing at boundaries. This function is used as a
     *  subroutine to inside_loop().
     */
    virtual int ray_cross( double x, double y ) const;

    bool geom_same( const MyDXFLine &line, double eps = 1.0e-6 ) const;

    virtual void scale( double s );
};


/*! \brief DXF arc entity class.
 *
 *
 */
class MyDXFArc : public MyDXFPathEntity
{

    MyDXFVec _pc;
    double _r;
    double _ang1; // Must be between 0 and 2 pi.
    double _ang2; // Must be between 0 and 2 pi.

    virtual void debug_print( std::ostream &os ) const;

public:

    MyDXFArc() : _r(1.0), _ang1(0.0), _ang2(2.0*M_PI) {};
    MyDXFArc( class MyDXFFile *dxf );
    virtual ~MyDXFArc() {}

    virtual MyDXFArc *copy( void ) const { return( new MyDXFArc( *this ) ); }

    virtual MyDXFVec start( void ) const {
	return( MyDXFVec(_pc[0] + _r*cos(_ang1), _pc[1] + _r*sin(_ang1), _pc[2] ) ); 
    }
    virtual MyDXFVec end( void ) const { 
	return( MyDXFVec(_pc[0] + _r*cos(_ang2), _pc[1] + _r*sin(_ang2), _pc[2] ) ); 
    }

    void set_pc( const MyDXFVec &pc ) { _pc = pc; }
    void set_r( double r ) { _r = r; }
    void set_ang1( double ang1 );
    void set_ang2( double ang2 );

    /*! \brief Reset arc according to end points.
     *
     *  Resets the center point and angles according to starting point
     *  \a s and ending point \a e.  Arc is assumed to go in
     *  counter-clockwise direction from start to end. Radius of arc
     *  is kept constant if distance from \a s to \a e is less than
     *  double the radius. Otherwise the center point is set to the
     *  middle of start and end and the radius is set to half of the
     *  distance of start to end.
     */
    void set_center_point( const MyDXFVec &s, const MyDXFVec &e );

    virtual void set_start( const MyDXFVec &s );
    virtual void set_end( const MyDXFVec &e );

    /*! \brief Check for ray crossing.
     *
     *  Check if ray going from point (x,y) downwards (negative y
     *  direction) crosses the entity. Return 1 if crosses odd number
     *  of times and 0 if even number of times. Return 2 in case of
     *  exact crossing at boundaries. This function is used as a
     *  subroutine to inside_loop().
     */
    virtual int ray_cross( double x, double y ) const;

    bool geom_same( const MyDXFArc &arc, double eps = 1.0e-6 ) const;

    virtual void scale( double s );
};


/*! \brief DXF circle entity class.
 *
 *
 */
class MyDXFCircle : public MyDXFPathEntity
{

    MyDXFVec _pc;
    double _r;

    virtual void debug_print( std::ostream &os ) const;

public:

    MyDXFCircle( class MyDXFFile *dxf );
    virtual ~MyDXFCircle() {}

    virtual MyDXFCircle *copy( void ) const { return( new MyDXFCircle( *this ) ); }

    virtual MyDXFVec start( void ) const { return( _pc+MyDXFVec(_r,0) ); }
    virtual MyDXFVec end( void ) const { return( _pc+MyDXFVec(_r,0) ); }

    virtual void set_start( const MyDXFVec &s ) {}
    virtual void set_end( const MyDXFVec &e ) {}

    /*! \brief Check for ray crossing.
     *
     *   Check if ray going from point (x,y) downwards (negative y
     *   direction) crosses the entity. Return 1 if crosses odd number
     *   of times and 0 if even number of times. Return 2 in case of
     *   exact crossing at boundaries. This function is used as a
     *   subroutine to inside_loop().
     */
    virtual int ray_cross( double x, double y ) const;

    bool geom_same( const MyDXFCircle &circle, double eps = 1.0e-6 ) const;

    virtual void scale( double s );
};


#define ATTACHMENT_POINT_TOP_LEFT       1
#define ATTACHMENT_POINT_TOP_CENTER     2
#define ATTACHMENT_POINT_TOP_RIGHT      3
#define ATTACHMENT_POINT_MIDDLE_LEFT    4
#define ATTACHMENT_POINT_MIDDLE_CENTER  5
#define ATTACHMENT_POINT_MIDDLE_RIGHT   6
#define ATTACHMENT_POINT_BOTTOM_LEFT    7
#define ATTACHMENT_POINT_BOTTOM_CENTER  8
#define ATTACHMENT_POINT_BOTTOM_RIGHT   9

#define DRAWING_DIRECTION_LEFT_TO_RIGHT 1
#define DRAWING_DIRECTION_TOP_TO_BOTTOM 3
#define DRAWING_DIRECTION_BY_STYLE      5


/*! \brief DXF text entity class.
 *
 *
 */
class MyDXFMText : public MyDXFEntity
{

    std::string _text;
    MyDXFVec    _p;
    double      _text_height;
    double      _rect_width;
    int16_t     _attachment_point;
    int16_t     _drawing_direction;

    virtual void debug_print( std::ostream &os ) const;

public:

    MyDXFMText( class MyDXFFile *dxf );
    virtual ~MyDXFMText() {}

    virtual MyDXFMText *copy( void ) const { return( new MyDXFMText( *this ) ); }

    virtual void scale( double s );
};



/*! \brief DXF entity selection.
 *
 *  A selection is a list of indexes, which point to entities in
 *  MyDXFEntities.
 */
class MyDXFEntitySelection
{

    std::vector<uint32_t> _selection;

public:

    MyDXFEntitySelection() {}
    ~MyDXFEntitySelection() {}

    uint32_t size() const { return( _selection.size() ); }
    void add_entity( uint32_t a ) { _selection.push_back( a ); }

    const uint32_t &operator()( int a ) const { return( _selection[a] ); }
    uint32_t &operator()( int a ) { return( _selection[a] ); }

    friend std::ostream &operator<<( std::ostream &os, const MyDXFEntitySelection &sel );
};


/*! \brief DXF entity database.
 *
 *  A database of entities.
 */
class MyDXFEntities
{

    std::vector<MyDXFEntity *> _entities;

public:


    MyDXFEntities() {}

    /*! \brief Construct new entities containing copies of selected intities in \a ent.
     */
    MyDXFEntities( MyDXFEntities *ent, MyDXFEntitySelection *sel );
    MyDXFEntities( class MyDXFFile *dxf );
    ~MyDXFEntities();

    /*
     *  ACCESS
     */

    uint32_t size() const { return( _entities.size() ); }
    const MyDXFEntity *get_entity( uint32_t a ) const { return( _entities[a] ); }
    MyDXFEntity *get_entity( uint32_t a ) { return( _entities[a] ); }

    /*
     *  GETTING SELECTION
     */

    /*! \brief Make a new selection containg all entities from database.
     */
    MyDXFEntitySelection *selection_all( void ) const;

    /*! \brief Make a new selection containg entities from named layer.
     */
    MyDXFEntitySelection *selection_layer( const std::string &layername ) const;

    /*! \brief Build complete loops.
     *
     * Make a subselection of a selection. The new selection will
     * contain only objects, which make up one or several complete
     * loops. Ending point of one entity is ensured to be exactly the
     * starting point of another entity. Errors of the size eps are
     * accepted and fixed between end-to-end matching.
     */
    MyDXFEntitySelection *selection_path_loop( MyDXFEntitySelection *selection,
					       double eps = 1.0e-6 );





    /*! \brief Check if two entities are geometrically same.
     *
     *   Checks if entity \a a is the geometrically same as entity \a
     *   b within error limit \a eps.
     */
    bool geom_same( uint32_t a, uint32_t b, double eps = 1.0e-6 ) const;

    /* ! \brief Check if point is inside a loop defined by a selection
     *   of entities.
     *
     *   The check is done assuming a 2D drawing in xy-plane. The
     *   check is done using ray shooting algorithm. If exact crossing
     *   happens perturbation algorithm is used. New test is performed
     *   at eps distance from the first.
     */
    bool inside_loop( MyDXFEntitySelection *selection, double x, double y, double eps = 1.0e-6 );





    void translate( MyDXFEntitySelection *selection, double dx, double dy, double dz );
    void rotate_x( MyDXFEntitySelection *selection, double y, double z, double ang );
    void rotate_y( MyDXFEntitySelection *selection, double x, double z, double ang );
    void rotate_z( MyDXFEntitySelection *selection, double x, double y, double ang );

    /* ! \brief Scale selected entities by factor s.
     */
    void scale( MyDXFEntitySelection *selection, double s );
};





#endif
