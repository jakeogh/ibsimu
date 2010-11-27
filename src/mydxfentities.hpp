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

#ifndef MY_DXF_ENTITIES_HPP
#define MY_DXF_ENTITIES_HPP 1


#include <stdint.h>
#include <vector>
#include <cairo.h>
#include "mydxffile.hpp"
#include "vec3d.hpp"
#include "transformation.hpp"


/*! \brief Entity type.
 */
enum EntityType {
    ENTITY_UNKNOWN = 0,
    ENTITY_LINE,
    ENTITY_LWPOLYLINE,
    ENTITY_ARC,
    ENTITY_CIRCLE,
    ENTITY_MTEXT,
    ENTITY_INSERT
};


/*! \brief DXF entity base class.
 *
 *  A general base class for all DXF entities. Contains data fields
 *  common to all entities.
 */
class MyDXFEntity
{

protected: 

    std::string _handle;
    std::string _layer;

    MyDXFEntity();

    //MyDXFEntity( const MyDXFEntity &ent );

    /*! \brief Propose a point to bounding box.
     *
     *  Updates bounding box value at min and max by including point p
     *  in the bounding box.
     */
    static void bbox_ppoint( Vec3D &min, Vec3D &max, const Vec3D &p );

    void write_common( class MyDXFFile *dxf, std::ofstream &ostr );
    void process_group( class MyDXFFile *dxf );
    void constructor_debug_print( void ) const;
    void debug_print_base( std::ostream &os ) const;
    virtual void debug_print( std::ostream &os ) const = 0;

public:

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFEntity() {}

    /*! \brief Get a new copy of entity.
     */
    virtual MyDXFEntity *copy( void ) const = 0;

    /*! \brief Explode into entities.
     *
     *  Break entity into atomic entities and tranform entities them
     *  with tranformation \a t. Add the tranformed entities to the
     *  database \a ent.
     */
    virtual void explode( class MyDXFEntities *ent, MyDXFFile *dxf, const Transformation *t ) const = 0;

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr ) = 0;

    /*! \brief Scale entity by factor \a s.
     */
    virtual void scale( class MyDXFFile *dxf, double s ) = 0;

    /*! \brief Set layer.
     */
    void set_layer( const std::string &layer ) { _layer = layer; }

    /*! \brief Get layer.
     */
    std::string get_layer( void ) const { return( _layer ); }

    /*! \brief Get entity type.
     */
    virtual EntityType get_type( void ) const = 0;

    /*! \brief Set entity handle.
     */
    void set_handle( const std::string &handle ) { _handle = handle; }

    /*! \brief Get entity handle.
     */
    std::string get_handle( void ) const { return( _handle ); }

    /*! \brief Plot entity with cairo
     *
     *  Plot the entity using the transformation \a from from the
     *  object space to cairo coordinates. The visible range is
     *  specified by \a range (xmin,ymin,xmax,ymax) in cairo
     *  coordinates.
     */
    virtual void plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		       const Transformation *t, const double range[4] ) const = 0;

    /*! \brief Return bounding box of entity
     */
    virtual void get_bbox( Vec3D &min, Vec3D &max, 
			   const class MyDXFFile *dxf, const Transformation *t ) const = 0;

    friend std::ostream &operator<<( std::ostream &os, const MyDXFEntity &ent );
};


/*! \brief DXF path entity base class.
 *
 *  A base class for two dimensional DXF entities, which can be part
 *  of a path. All path entities have a start point and an end point,
 *  that can be read and set.
 */
class MyDXFPathEntity : public MyDXFEntity
{

protected: 

    MyDXFPathEntity() {}

    MyDXFPathEntity( const MyDXFEntity &ent ) : MyDXFEntity(ent) {}

public:

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFPathEntity() {}

    /*! \brief Get start point of path entity.
     */
    virtual Vec3D start( void ) const = 0;

    /*! \brief Get end point of path entity.
     */
    virtual Vec3D end( void ) const = 0;

    /*! \brief Set start point of path entity.
     */
    virtual void set_start( const Vec3D &s ) = 0;

    /*! \brief Set end point of path entity.
     */
    virtual void set_end( const Vec3D &e ) = 0;

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

    Vec3D _p1;
    Vec3D _p2;

    virtual void debug_print( std::ostream &os ) const;

public:

    /*! \brief Default constructor.
     */
    MyDXFLine() {}

    /*! \brief Constructor for copying MyDXFEntity properties.
     */
    MyDXFLine( const MyDXFEntity &ent ) : MyDXFPathEntity(ent) {}

    /*! \brief Construct line entity by reading from DXF file.
     */
    MyDXFLine( class MyDXFFile *dxf );

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFLine() {}

    /*! \brief Get a new copy of entity.
     */
    virtual MyDXFLine *copy( void ) const { return( new MyDXFLine( *this ) ); }

    /*! \brief Explode into entities.
     *
     *  Break entity into atomic entities and tranform entities them
     *  with tranformation \a t. Add the tranformed entities to the
     *  database \a ent.
     */
    virtual void explode( class MyDXFEntities *ent, MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Get entity type.
     */
    virtual EntityType get_type( void ) const { return( ENTITY_LINE ); }

    /*! \brief Get start point of path entity.
     */
    virtual Vec3D start( void ) const { return( _p1 ); }

    /*! \brief Get end point of path entity.
     */
    virtual Vec3D end( void ) const { return( _p2 ); }

    /*! \brief Set start point of path entity.
     */
    virtual void set_start( const Vec3D &s );

    /*! \brief Set end point of path entity.
     */
    virtual void set_end( const Vec3D &e );

    /*! \brief Check for ray crossing.
     *
     *  Check if ray going from point (x,y) downwards (negative y
     *  direction) crosses the entity. Return 1 if crosses odd number
     *  of times and 0 if even number of times. Return 2 in case of
     *  exact crossing at boundaries. This function is used as a
     *  subroutine to inside_loop().
     */
    virtual int ray_cross( double x, double y ) const;

    /*! \brief Check if two entities are geometrically same.
     *
     *   Checks if entity \a a is the geometrically same as entity \a
     *   b within error limit \a eps.
     */
    bool geom_same( const MyDXFLine &line, double eps = 1.0e-6 ) const;

    /*! \brief Plot entity with cairo
     *
     *  Plot the entity using the transformation \a from from the
     *  object space to cairo coordinates. The visible range is
     *  specified by \a range (xmin,ymin,xmax,ymax) in cairo
     *  coordinates.
     */
    virtual void plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		       const Transformation *t, const double range[4] ) const;

   /*! \brief Return bounding box of entity
     */
    virtual void get_bbox( Vec3D &min, Vec3D &max, 
			   const class MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Scale entity by factor \a s.
     */
    virtual void scale( class MyDXFFile *dxf, double s );
};


/* Polyline flags */
#define LWPOLYLINE_CLOSED_MASK  1
#define LWPOLYLINE_LINEGEN_MASK 128


/*! \brief DXF LW polyline entity class.
 *
 *  
 */
class MyDXFLWPolyline : public MyDXFPathEntity
{

    std::vector<Vec3D> _p;       /*!< \brief Vector with x and y coordinates of 
				  *   vertex and bulge as z coordinate */
    int16_t            _flags;

    virtual void debug_print( std::ostream &os ) const;

public:

    /*! \brief Construct line entity by reading from DXF file.
     */
    MyDXFLWPolyline( class MyDXFFile *dxf );

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFLWPolyline() {}

    /*! \brief Get a new copy of entity.
     */
    virtual MyDXFLWPolyline *copy( void ) const { return( new MyDXFLWPolyline( *this ) ); }

    /*! \brief Explode into entities.
     *
     *  Break entity into atomic entities and tranform entities them
     *  with tranformation \a t. Add the tranformed entities to the
     *  database \a ent.
     */
    virtual void explode( class MyDXFEntities *ent, MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Get entity type.
     */
    virtual EntityType get_type( void ) const { return( ENTITY_LWPOLYLINE ); }

    /*! \brief Get start point of path entity.
     */
    virtual Vec3D start( void ) const { return( _p[0] ); }

    /*! \brief Get end point of path entity.
     */
    virtual Vec3D end( void ) const { return( _p[_p.size()-1] ); }

    /*! \brief Set start point of path entity.
     */
    virtual void set_start( const Vec3D &s );

    /*! \brief Set end point of path entity.
     */
    virtual void set_end( const Vec3D &e );

    /*! \brief Check for ray crossing.
     *
     *  Check if ray going from point (x,y) downwards (negative y
     *  direction) crosses the entity. Return 1 if crosses odd number
     *  of times and 0 if even number of times. Return 2 in case of
     *  exact crossing at boundaries. This function is used as a
     *  subroutine to inside_loop().
     */
    virtual int ray_cross( double x, double y ) const;

    /*! \brief Check if two entities are geometrically same.
     *
     *   Checks if entity \a a is the geometrically same as entity \a
     *   b within error limit \a eps.
     */
    bool geom_same( const MyDXFLWPolyline &line, double eps = 1.0e-6 ) const;

    /*! \brief Get number of vertices in entity.
     */
    uint32_t size() const { return( _p.size() ); }

    /*! \brief Get vertix \a i.
     */
    Vec3D vertex( uint32_t i ) const { return( _p[i] ); }

    /*! \brief Is entity closed path?
     */
    bool closed( void ) const { return( _flags & LWPOLYLINE_CLOSED_MASK ); }

    /*! \brief Plot entity with cairo
     *
     *  Plot the entity using the transformation \a from from the
     *  object space to cairo coordinates. The visible range is
     *  specified by \a range (xmin,ymin,xmax,ymax) in cairo
     *  coordinates.
     */
    virtual void plot( const class MyDXFFile *dxf, cairo_t *cairo,
		       const Transformation *t, const double range[4] ) const;

    /*! \brief Return bounding box of entity.
     */
    virtual void get_bbox( Vec3D &min, Vec3D &max, 
			   const class MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Scale entity by factor \a s.
     */
    virtual void scale( class MyDXFFile *dxf, double s );
};


/*! \brief DXF arc entity class.
 *
 *  An arc entity is defined with a center point of the arc, radius,
 *  starting angle and ending angle.
 */
class MyDXFArc : public MyDXFPathEntity
{

    Vec3D    _pc;
    double   _r;
    double   _ang1; // Must be between 0 and 2 pi.
    double   _ang2; // Must be between 0 and 2 pi.

    virtual void debug_print( std::ostream &os ) const;

public:

    /*! \brief Default constructor.
     */
    MyDXFArc() : _r(1.0), _ang1(0.0), _ang2(2.0*M_PI) {};

    /*! \brief Construct arc entity by reading from DXF file.
     */
    MyDXFArc( class MyDXFFile *dxf );

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFArc() {}

    /*! \brief Get a new copy of entity.
     */
    virtual MyDXFArc *copy( void ) const { return( new MyDXFArc( *this ) ); }

    /*! \brief Explode into entities.
     *
     *  Break entity into atomic entities and tranform entities them
     *  with tranformation \a t. Add the tranformed entities to the
     *  database \a ent.
     */
    virtual void explode( class MyDXFEntities *ent, MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Get entity type.
     */
    virtual EntityType get_type( void ) const { return( ENTITY_ARC ); }

    /*! \brief Get center point of arc.
     */
    Vec3D center( void ) const { return( _pc ); }

    /*! \brief Get radius of arc.
     */
    double radius( void ) const { return( _r ); }

    /*! \brief Get start point of path entity.
     */
    virtual Vec3D start( void ) const {
	return( Vec3D(_pc[0] + _r*cos(_ang1), _pc[1] + _r*sin(_ang1), _pc[2] ) ); 
    }

    /*! \brief Get end point of path entity.
     */
    virtual Vec3D end( void ) const { 
	return( Vec3D(_pc[0] + _r*cos(_ang2), _pc[1] + _r*sin(_ang2), _pc[2] ) ); 
    }

    /*! \brief Set center point.
     */
    void set_pc( const Vec3D &pc ) { _pc = pc; }

    /*! \brief Set radius.
     */
    void set_r( double r ) { _r = r; }

    /*! \brief Set start angle.
     */
    void set_ang1( double ang1 );

    /*! \brief Set end angle.
     */
    void set_ang2( double ang2 );
    
    /*! \brief Get start angle.
     */
    double get_ang1( void ) const { return( _ang1 ); }

    /*! \brief Get end angle.
     */
    double get_ang2( void ) const { return( _ang2 ); }

    /*! \brief Set arc according to center and end points.
     *
     *  End point is only used for the angle.
     */
    void set_center_and_ends( const Vec3D &c, const Vec3D &s, const Vec3D &e );

    /*! \brief Set arc according to end points.
     *
     *  Resets the center point and angles according to starting point
     *  \a s and ending point \a e.  Arc is assumed to go in
     *  counter-clockwise direction from start to end. Radius of arc
     *  is kept constant if distance from \a s to \a e is less than
     *  double the radius. Otherwise the center point is set to the
     *  middle of start and end and the radius is set to half of the
     *  distance of start to end.
     */
    void set_center_point( const Vec3D &s, const Vec3D &e );

    /*! \brief Set start point of path entity.
     */
    virtual void set_start( const Vec3D &s );

    /*! \brief Set end point of path entity.
     */
    virtual void set_end( const Vec3D &e );

    /*! \brief Check for ray crossing.
     *
     *  Check if ray going from point (x,y) downwards (negative y
     *  direction) crosses the entity. Return 1 if crosses odd number
     *  of times and 0 if even number of times. Return 2 in case of
     *  exact crossing at boundaries. This function is used as a
     *  subroutine to inside_loop().
     */
    virtual int ray_cross( double x, double y ) const;

    /*! \brief Check if two entities are geometrically same.
     *
     *   Checks if entity \a a is the geometrically same as entity \a
     *   b within error limit \a eps.
     */
    bool geom_same( const MyDXFArc &arc, double eps = 1.0e-6 ) const;

    /*! \brief Plot entity with cairo
     *
     *  Plot the entity using the transformation \a from from the
     *  object space to cairo coordinates. The visible range is
     *  specified by \a range (xmin,ymin,xmax,ymax) in cairo
     *  coordinates.
     */
    virtual void plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		       const Transformation *t, const double range[4] ) const;

    /*! \brief Return bounding box of entity
     */
    virtual void get_bbox( Vec3D &min, Vec3D &max, 
			   const class MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Scale entity by factor \a s.
     */
    virtual void scale( class MyDXFFile *dxf, double s );
};


/*! \brief DXF circle entity class.
 *
 *
 */
class MyDXFCircle : public MyDXFPathEntity
{

    Vec3D    _pc;
    double   _r;

    virtual void debug_print( std::ostream &os ) const;

public:

    /*! \brief Default constructor.
     */
    MyDXFCircle() : _r(1.0) {};

    /*! \brief Construct circle entity by reading from DXF file.
     */
    MyDXFCircle( class MyDXFFile *dxf );

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFCircle() {}

    /*! \brief Get a new copy of entity.
     */
    virtual MyDXFCircle *copy( void ) const { return( new MyDXFCircle( *this ) ); }

    /*! \brief Explode into entities.
     *
     *  Break entity into atomic entities and tranform entities them
     *  with tranformation \a t. Add the tranformed entities to the
     *  database \a ent.
     */
    virtual void explode( class MyDXFEntities *ent, MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Get entity type.
     */
    virtual EntityType get_type( void ) const { return( ENTITY_CIRCLE ); }

    /*! \brief Get center point of circle.
     */
    Vec3D center( void ) const { return( _pc ); }

    /*! \brief Get radius of circle.
     */
    double radius( void ) const { return( _r ); }

    /*! \brief Get start point of path entity.
     */
    virtual Vec3D start( void ) const { return( _pc+Vec3D(_r,0) ); }

    /*! \brief Get end point of path entity.
     */
    virtual Vec3D end( void ) const { return( _pc+Vec3D(_r,0) ); }

    /*! \brief Set start point of path entity.
     */
    virtual void set_start( const Vec3D &s ) {}

    /*! \brief Set end point of path entity.
     */
    virtual void set_end( const Vec3D &e ) {}

    /*! \brief Check for ray crossing.
     *
     *   Check if ray going from point (x,y) downwards (negative y
     *   direction) crosses the entity. Return 1 if crosses odd number
     *   of times and 0 if even number of times. Return 2 in case of
     *   exact crossing at boundaries. This function is used as a
     *   subroutine to inside_loop().
     */
    virtual int ray_cross( double x, double y ) const;

    /*! \brief Check if two entities are geometrically same.
     *
     *   Checks if entity \a a is the geometrically same as entity \a
     *   b within error limit \a eps.
     */
    bool geom_same( const MyDXFCircle &circle, double eps = 1.0e-6 ) const;

    /*! \brief Plot entity with cairo
     *
     *  Plot the entity using the transformation \a from from the
     *  object space to cairo coordinates. The visible range is
     *  specified by \a range (xmin,ymin,xmax,ymax) in cairo
     *  coordinates.
     */
    virtual void plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		       const Transformation *t, const double range[4] ) const;

    /*! \brief Return bounding box of entity
     */
    virtual void get_bbox( Vec3D &min, Vec3D &max, 
			   const class MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Scale entity by factor \a s.
     */
    virtual void scale( class MyDXFFile *dxf, double s );
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
    Vec3D       _p;
    double      _text_height;
    double      _rect_width;
    int16_t     _attachment_point;
    int16_t     _drawing_direction;

    virtual void debug_print( std::ostream &os ) const;

public:

    /*! \brief Default constructor.
     */
    MyDXFMText() : _text_height(1.0), _rect_width(1.0), 
		   _attachment_point(ATTACHMENT_POINT_TOP_LEFT),
		   _drawing_direction(DRAWING_DIRECTION_LEFT_TO_RIGHT) {};

    /*! \brief Construct line entity by reading from DXF file.
     */
    MyDXFMText( class MyDXFFile *dxf );

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFMText() {}

    /*! \brief Get a new copy of entity.
     */
    virtual MyDXFMText *copy( void ) const { return( new MyDXFMText( *this ) ); }

    /*! \brief Explode into entities.
     *
     *  Break entity into atomic entities and tranform entities them
     *  with tranformation \a t. Add the tranformed entities to the
     *  database \a ent.
     */
    virtual void explode( class MyDXFEntities *ent, MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Get entity type.
     */
    virtual EntityType get_type( void ) const { return( ENTITY_MTEXT ); }

    /*! \brief Plot entity with cairo
     *
     *  Plot the entity using the transformation \a from from the
     *  object space to cairo coordinates. The visible range is
     *  specified by \a range (xmin,ymin,xmax,ymax) in cairo
     *  coordinates.
     */
    virtual void plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		       const Transformation *t, const double range[4] ) const;

    /*! \brief Return bounding box of entity
     */
    virtual void get_bbox( Vec3D &min, Vec3D &max, 
			   const class MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Scale entity by factor \a s.
     */
    virtual void scale( class MyDXFFile *dxf, double s );
};



/*! \brief DXF insert entity class.
 *
 *
 */
class MyDXFInsert : public MyDXFEntity
{

    std::string _block_name;

    Vec3D       _p;        // Insertion point
    Vec3D       _scale;

    double      _rotation;

    int16_t     _col_count;
    int16_t     _row_count;

    double      _col_spacing;
    double      _row_spacing;

    virtual void debug_print( std::ostream &os ) const;

public:

    /*! \brief Default constructor.
     */
    MyDXFInsert();

    /*! \brief Construct entity by reading from DXF file.
     */
    MyDXFInsert( class MyDXFFile *dxf );

    /*! \brief Virtual destructor.
     */
    virtual ~MyDXFInsert() {}

    /*! \brief Get a new copy of entity.
     */
    virtual MyDXFInsert *copy( void ) const { return( new MyDXFInsert( *this ) ); }

    /*! \brief Explode into entities.
     *
     *  Break entity into atomic entities and tranform entities them
     *  with tranformation \a t. Add the tranformed entities to the
     *  database \a ent.
     */
    virtual void explode( class MyDXFEntities *ent, MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Write dxf file to stream.
     */
    virtual void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Get entity type.
     */
    virtual EntityType get_type( void ) const { return( ENTITY_INSERT ); }

    /*! \brief Plot entity with cairo
     *
     *  Plot the entity using the transformation \a from from the
     *  object space to cairo coordinates. The visible range is
     *  specified by \a range (xmin,ymin,xmax,ymax) in cairo
     *  coordinates.
     */
    virtual void plot( const class MyDXFFile *dxf, cairo_t *cairo, 
		       const Transformation *t, const double range[4] ) const;

    /*! \brief Return bounding box of entity
     */
    virtual void get_bbox( Vec3D &min, Vec3D &max, 
			   const class MyDXFFile *dxf, const Transformation *t ) const;

    /*! \brief Scale entity by factor \a s.
     */
    virtual void scale( class MyDXFFile *dxf, double s );

   /*! \brief Return name of the block inserted.
     */
    const std::string &block_name( void ) const { return( _block_name ); }
};



/*! \brief DXF entity selection.
 *
 *  %MyDXFEntitySelection object is a list of indexes, which point to
 *  entities in the MyDXFEntities.
 */
class MyDXFEntitySelection
{

    std::vector<uint32_t> _selection;

public:

    /*! \brief Construct empty selection.
     */
    MyDXFEntitySelection() {}

    /*! \brief Destructor.
     */
    ~MyDXFEntitySelection() {}

    /*! \brief Return number of entities in selection.
     */
    uint32_t size() const { return( _selection.size() ); }

    /*! \brief Add entity number \a a in selection.
     */
    void add_entity( uint32_t a ) { _selection.push_back( a ); }

    /*! \brief Get a const reference to entity number in selection at location \a a.
     */
    const uint32_t &operator()( int a ) const { 
	if( a < 0 || a >= (int)_selection.size() )
	    throw( Error( ERROR_LOCATION, "index out of range" ) );
	return( _selection[a] ); 
    }

    /*! \brief Get reference to entity number in selection at location \a a.
     */
    uint32_t &operator()( int a ) { 
	if( a < 0 || a >= (int)_selection.size() )
	    throw( Error( ERROR_LOCATION, "index out of range" ) );
	return( _selection[a] ); 
    }

    friend std::ostream &operator<<( std::ostream &os, const MyDXFEntitySelection &sel );
};


/*! \brief DXF entity database.
 *
 *  A database of entities. Also responsible for reading entities from
 *  a DXF file. All supported entities are saved to database. All
 *  others are silently ignored.
 */
class MyDXFEntities
{

    std::vector<MyDXFEntity *> _entities;

public:


    /*! \brief Construct empty entities database.
     */
    MyDXFEntities() {}

    /*! \brief Construct new entities containing copies of selected intities in \a ent.
     */
    MyDXFEntities( MyDXFEntities *ent, MyDXFEntitySelection *sel );

    /*! \brief Construct entities database by reading from DXF file.
     *
     *  Called with reading_block = true if called from inside BLOCKS
     *  section.
     */
    MyDXFEntities( class MyDXFFile *dxf, bool reading_blocks = false );

    /*! \brief Destructor.
     */
    ~MyDXFEntities();


    /*! \brief Write entities section of dxf file to stream.
     */
    void write( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Write a list of entities to stream.
     *
     *  This function only writes the entities withing the object. It
     *  can be called for writing the entities section or a block.
     */
    void write_entities( class MyDXFFile *dxf, std::ofstream &ostr );

    /*! \brief Return number of entities.
     */
    uint32_t size() const { return( _entities.size() ); }

    /*! \brief Return const pointer to entity \a a.
     */
    const MyDXFEntity *get_entity( uint32_t a ) const { return( _entities[a] ); }

    /*! \brief Return pointer to entity \a a.
     */
    MyDXFEntity *get_entity( uint32_t a ) { return( _entities[a] ); }



    
    /*! \brief Add entity to list.
     *
     *  No copy of entity is made. The pointer is saved to the database.
     */
    void add_entity( MyDXFEntity *e ) { _entities.push_back( e ); }


    /*! \brief Make a new selection containg all entities from database.
     */
    MyDXFEntitySelection *selection_all( void ) const;

    /*! \brief Make a new selection containg entities from named layer.
     */
    MyDXFEntitySelection *selection_layer( const std::string &layername ) const;

    /*! \brief Make a new selection containg entities of given type.
     */
    MyDXFEntitySelection *selection_type( EntityType type ) const;

    /*! \brief Build complete loops.
     *
     * Make a new subselection of a selection. The new selection will
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

    /*! \brief Plot selected entities with cairo
     *
     *  Plot the entities using the transformation \a t from the
     *  object space to cairo coordinates. The visible range is
     *  specified by \a range (xmin,ymin,xmax,ymax) in cairo
     *  coordinates.
     *
     *  Selection can be a NULL pointer to plot all entities.
     */
    void plot( const MyDXFEntitySelection *selection, const class MyDXFFile *dxf, 
	       cairo_t *cairo, const Transformation *t, const double range[4] ) const;

    /*! \brief Get bounding box containing all entities in selection.
     */
    void get_bbox( const MyDXFEntitySelection *selection, Vec3D &min, Vec3D &max, 
		   const class MyDXFFile *dxf, const Transformation *t ) const;



    /*
    void translate( MyDXFEntitySelection *selection, double dx, double dy, double dz );
    void rotate_x( MyDXFEntitySelection *selection, double y, double z, double ang );
    void rotate_y( MyDXFEntitySelection *selection, double x, double z, double ang );
    void rotate_z( MyDXFEntitySelection *selection, double x, double y, double ang );
    */

    /*! \brief Scale selected entities by factor s.
     *
     *  Selection can be a NULL pointer to plot all entities.
     */
    void scale( MyDXFEntitySelection *selection, class MyDXFFile *dxf, double s );


    /*! \brief Remove selected entities.
     *
     *  Selection can be a NULL pointer to remove all entities. The
     *  selection is invalid after this operation and should not be
     *  used further. Also all other selections are invalidated by
     *  this operation because entity indices change.
     */
    void remove( MyDXFEntitySelection *selection );


    /*! \brief Explode selected insert entities.
     *
     *  The insert entities are expoded to contain just primitive
     *  entities with no dependencies to blocks. Selection can be a
     *  NULL pointer to explode all entities.
     */
    void explode( MyDXFEntitySelection *selection, class MyDXFFile *dxf );

    /*! \brief Explode all entities to \a ent.
     *
     *  Explode and add all entities into ent using transformation t.
     */
    void explode( MyDXFEntities *ent, class MyDXFFile *dxf, const Transformation *t ) const;


    /*! \brief Print debugging information to os.
     */
    void debug_print( std::ostream &os ) const;


};





#endif



