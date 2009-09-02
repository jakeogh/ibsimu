#ifndef FIELDDIAGPLOT_HPP
#define FIELDDIAGPLOT_HPP 1


#include "frame.hpp"
#include "xygraph.hpp"
#include "vec3d.hpp"
#include "geometry.hpp"
#include "scalarfield.hpp"
#include "efield.hpp"
#include "vectorfield.hpp"


enum field_diag_type_e {
    FIELDD_DIAG_NONE = 0,
    FIELDD_DIAG_EPOT,
    FIELDD_DIAG_EFIELD,
    FIELDD_DIAG_EFIELD_X,
    FIELDD_DIAG_EFIELD_Y,
    FIELDD_DIAG_EFIELD_Z,
    FIELDD_DIAG_SCHARGE,
    FIELDD_DIAG_BFIELD,
    FIELDD_DIAG_BFIELD_X,
    FIELDD_DIAG_BFIELD_Y,
    FIELDD_DIAG_BFIELD_Z
};


enum field_loc_type_e {
    FIELDD_LOC_NONE = 0,
    FIELDD_LOC_X,
    FIELDD_LOC_Y,
    FIELDD_LOC_Z,
    FIELDD_LOC_DIST
};


class FieldDiagPlot {

    Frame              *_frame;

    const Geometry     *_geom;
    const ScalarField  *_epot;
    const Efield       *_efield;
    const ScalarField  *_scharge;
    const VectorField  *_bfield;

    size_t              _N;
    Vec3D               _x1;
    Vec3D               _x2;

    field_diag_type_e   _diag[2];
    field_loc_type_e    _loc[2];

    XYGraph            *_graph[2];

public:

    FieldDiagPlot( Frame *frame, const Geometry *geom );

    ~FieldDiagPlot();

    void set_epot( const ScalarField *epot ) {
	_epot = epot;
    }

    void set_efield( const Efield *efield ) {
	_efield = efield;
    }

    void set_scharge( const ScalarField *scharge ) {
	_scharge = scharge;
    }

    void set_bfield( const VectorField *bfield ) {
	_bfield = bfield;
    }

    void set_coordinates( size_t N, const Vec3D &x1, const Vec3D &x2 ) {
	_N = N;
	_x1 = x1;
	_x2 = x2;
    }

    void set_diagnostic( const field_diag_type_e diag[2], const field_loc_type_e loc[2] ) {
	_diag[0] = diag[0];
	_diag[1] = diag[1];
	_loc[0] = loc[0];
	_loc[1] = loc[1];
    }

    /*! \brief Rebuild plot.
     */
    void build_plot( void );

};



#endif
