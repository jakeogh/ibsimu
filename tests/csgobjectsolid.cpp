/*! \file csgobjectsolid.cpp 
 *  \test Test for CSGObjectSolid.
 */


#include <iostream>
#include <iomanip>
#include "bicgstab_solver.hpp"
#include "epot_problem.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "csgobject_solid.hpp"
#include "epot_efield.hpp"
#include "geomplotter.hpp"
#include "error.hpp"


using namespace std;


CSGObject *object1( void )
{
    CSGSphere *sphere1 = new CSGSphere();
    sphere1->translate( 0.7, 0, 0 );

    CSGSphere *sphere2 = new CSGSphere();
    //sphere2->translate( 0, 0, 2 );

    CSGDifference *diff = new CSGDifference( );
    diff->add( sphere1 );
    diff->add( sphere2 );

    return( diff );
}


void test( void )
{
    Geometry geom( MODE_2D, Int3D(101,101,1), Vec3D(-2.5,-2.5,0), 0.05 );
    CSGObject *o1 = object1();
    Solid *s1 = new CSGObjectSolid( o1 );
    geom.set_solid( 7, s1 );
    geom.set_boundary( 1,  Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 2,  Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 3,  Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4,  Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7,  Bound(BOUND_DIRICHLET,  0.0) );
    geom.build_mesh();

    GeomPlotter gplotter( &geom );
    gplotter.set_mesh( true );
    gplotter.plot_png( "csgobjectsolid.png" );
}


int main( void )
{
    try {
        test();
    } catch ( Error e ) {
        cout << "Error in " << e._loc._file << ":" << e._loc._line 
             << " in " << e._loc._func << "(): " << e._error_str << "\n";
        exit( 1 );
    }

    return( 0 );
}


