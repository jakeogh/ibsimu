/*! \file particles.hpp
 *  \brief Header file defining particle point and particle classes.
 */


#ifndef PARTICLES_HPP
#define PARTICLES_HPP 1


#include <vector>
#include <gsl/gsl_errno.h>
#include "geometry.hpp"
#include "scalarfield.hpp"
#include "vectorfield.hpp"
#include "efield.hpp"
#include "vec3d.hpp"


/* Integer error value that is supposed to diffed from internal GSL error values */
#define IBSIMU_DERIV_ERROR 201


/*! \brief %Particle status enum.
 *
 *  %Particle status can be either OK (\a PARTICLE_OK), out of geometry
 *  (\a PARTICLE_OUT), collided (\a PARTICLE_COLL), badly defined (\a
 *  PARTICLE_BADDEF), simulated to time limit (\a PARTICLE_TIME) or
 *  simulated to step limit (\a PARTICLE_NSTP).
 */
enum particle_status_e {
    PARTICLE_OK = 0,
    PARTICLE_OUT,
    PARTICLE_COLL,
    PARTICLE_BADDEF,
    PARTICLE_TIME,
    PARTICLE_NSTP 
};


/* Atomic mass unit */
#define MASS_U        1.66053873e-27
/* Elementary charge */
#define CHARGE_E      1.602176462e-19




/* ************************************************************************************************ *
 * Particle point classes                                                                           *
 * ************************************************************************************************ */


/*! \brief %Particle point base class
 *
 *  A parent class for a common handle to %Particle point classes.
 */
class ParticlePBase
{

};


/*! \brief %Particle point class for 2D.
 *
 *  A collection of particle coordinates (t,x,vx,y,vy) to describe
 *  time, location and velocity in 2D.
 */
class ParticleP2D : public ParticlePBase 
{
    double         _x[5];   /*!< \brief Point coordinates */

public:

    /*! \brief Default constuctor.
     */
    ParticleP2D() {}

    /*! \brief Constructor for 2D particle point.
     */
    ParticleP2D( double t, double x, double vx, double y, double vy ) {
	_x[0] = t; _x[1] = x; _x[2] = vx; _x[3] = y; _x[4] = vy;
    }

    /*! \brief Returns geometry mode.
     */
    static geom_mode_e geom_mode() { return(MODE_2D); }

    /*! \brief Returns number of dimensions for geometry.
     */
    static size_t dim() { return(2); }

    /*! \brief Returns number of coordinates used for particle point.
     */
    static size_t size() { return(5); }

    /*! \brief Returns time derivatives \a dxdt of coordinates at time
     *  \a t and coordinates \a x = (x,vx,y,vy) for one particle.
     *
     *  The calculation of particle trajectory is done by integrating
     *  the Lorentz equation in a form of a set of ordinary
     *  differential equations. In the case of 2D the set is:
     *  \f{eqnarray*} \frac{dx}{dt}   &=& v_x \\ 
                      \frac{dy}{dt}   &=& v_y \\ 
                      \frac{dv_x}{dt} &=& a_x = \frac{q}{m}(E_x + v_y B_z) \\ 
                      \frac{dv_y}{dt} &=& a_y = \frac{q}{m}(E_y - v_x B_z) \f}
     */
    static int get_derivatives( double t, const double *x, double *dxdt, void *data );

    /*! \brief Return the number of trajectory intersections with
     *  plane \a crd = val on the trajectory from \a x1 to \a
     *  x2. Intersection points are appended to vector \a intsc.
     */
    static int trajectory_intersections_at_plane( std::vector<ParticleP2D> &intsc, 
						  int crd, double val,
						  const ParticleP2D &x1, const ParticleP2D &x2 );

    /*! \brief Returns the location of particle point in Vec3D.
     */
    Vec3D location() const { return( Vec3D( _x[1], _x[3], 0.0 ) ); }

    /*! \brief Returns the velocity of particle point in Vec3D.
     */
    Vec3D velocity() const { return( Vec3D( _x[2], _x[4], 0.0 ) ); }

    /*! \brief Returns speed of particle.
     */
    double speed() { return( sqrt(_x[2]*_x[2] + _x[4]*_x[4]) ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator[]( int i ) { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator[]( int i ) const { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator()( int i ) { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator()( int i ) const { return( _x[i] ); }

    ParticleP2D operator+( const ParticleP2D &pp ) const { 
	ParticleP2D res;
	res[0] = _x[0] + pp[0];
	res[1] = _x[1] + pp[1];
	res[2] = _x[2] + pp[2];
	res[3] = _x[3] + pp[3];
	res[4] = _x[4] + pp[4];
	return( res );
    }

    ParticleP2D operator-( const ParticleP2D &pp ) const { 
	ParticleP2D res;
	res[0] = _x[0] - pp[0];
	res[1] = _x[1] - pp[1];
	res[2] = _x[2] - pp[2];
	res[3] = _x[3] - pp[3];
	res[4] = _x[4] - pp[4];
	return( res );
    }

    ParticleP2D operator*( double x ) const { 
	ParticleP2D res;
	res[0] = _x[0]*x;
	res[1] = _x[1]*x;
	res[2] = _x[2]*x;
	res[3] = _x[3]*x;
	res[4] = _x[4]*x;
	return( res );
    }

    friend ParticleP2D operator*( double x, const ParticleP2D &pp );
};


inline std::ostream &operator<<( std::ostream &os, const ParticleP2D &pp )
{
    os << std::setw(12) << pp(0) << " "
       << std::setw(12) << pp(1) << " "
       << std::setw(12) << pp(2) << " "
       << std::setw(12) << pp(3) << " "
       << std::setw(12) << pp(4);
    return( os );
}


inline ParticleP2D operator*( double x, const ParticleP2D &pp )
{
    ParticleP2D res;
    res[0] = pp[0]*x;
    res[1] = pp[1]*x;
    res[2] = pp[2]*x;
    res[3] = pp[3]*x;
    res[4] = pp[4]*x;
    return( res );
}


/*! \brief %Particle point class for cylindrical coordinates.
 *
 *  A collection of particle coordinates (t,x,vx,y,vy,w) to
 *  describe time, location and velocity in cylindrical coordinates.
 */
class ParticlePCyl : public ParticlePBase 
{
    double         _x[6];   /*!< \brief Point coordinates */

public:

    /*! \brief Default constuctor.
     */
    ParticlePCyl() {}

    /*! \brief Constructor for cylindrical particle point.
     */
    ParticlePCyl( double t, double x, double vx, double r, double vr, double w ) {
	_x[0] = t; _x[1] = x; _x[2] = vx; _x[3] = r; _x[4] = vr; _x[5] = w;
    }

    /*! \brief Returns geometry mode.
     */
    static geom_mode_e geom_mode() { return(MODE_CYL); }

    /*! \brief Returns number of dimensions for geometry.
     */
    static size_t dim() { return(2); }

    /*! \brief Returns number of coordinates used for particle point.
     */
    static size_t size() { return(6); }

    /*! \brief Returns time derivatives \a dxdt of coordinates at time
     *  \a t and coordinates \a x = (x,vx,r,vr,w) for one particle.
     *
     *  The calculation of particle trajectory is done by integrating
     *  the Lorentz equation in a form of a set of ordinary
     *  differential equations. In the case of cylindrical coordinates 
     *  the set is:
     *  \f{eqnarray*} \frac{dx}{dt}   &=& v_x \\ 
                      \frac{dr}{dt}   &=& v_r \\ 
                      \frac{dv_x}{dt} &=& a_x = \frac{q}{m}(E_x + v_r B_{\theta} - v_\theta B_r) \\ 
                      \frac{dv_r}{dt} &=& a_r + r\left(\frac{d\theta}{dt}\right)^2 = 
                        \frac{q}{m}(E_y + v_{\theta} B_x - v_x B_{\theta}) + r \left(\frac{d\theta}{dt}\right)^2 \\
                      \frac{d^2\theta}{dt^2} &=& \frac{1}{r} \left(a_\theta-\frac{dr}{dt} \frac{d\theta}{dt}\right) =
		        \frac{1}{r} \left(\frac{q}{m}(v_x B_r - v_r B_x)-2 \frac{dr}{dt} \frac{d\theta}{dt}\right), \f}
     *  where \f$ v_\theta = r \frac{d\theta}{dt} \f$.
     */
    static int get_derivatives( double t, const double *x, double *dxdt, void *data );

    /*! \brief Return the number of trajectory intersections with
     *  plane \a crd = val on the trajectory from \a x1 to \a
     *  x2. Intersection points are appended to vector \a intsc.
     */
    static int trajectory_intersections_at_plane( std::vector<ParticlePCyl> &intsc, 
						  int crd, double val,
						  const ParticlePCyl &x1, const ParticlePCyl &x2 );

    /*! \brief Returns the location of particle point in Vec3D.
     */
    Vec3D location() const { return( Vec3D( _x[1], _x[3], 0.0 ) ); }

    /*! \brief Returns the velocity of particle point in Vec3D.
     */
    Vec3D velocity() const { return( Vec3D( _x[2], _x[4], _x[5]*_x[3] ) ); }

    /*! \brief Returns speed of particle.
     */
    double speed() { return( sqrt(_x[2]*_x[2] + _x[4]*_x[4] + _x[3]*_x[3]*_x[5]*_x[5]) ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator[]( int i ) { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator[]( int i ) const { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator()( int i ) { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator()( int i ) const { return( _x[i] ); }

    ParticlePCyl operator+( const ParticlePCyl &pp ) const { 
	ParticlePCyl res;
	res[0] = _x[0] + pp[0];
	res[1] = _x[1] + pp[1];
	res[2] = _x[2] + pp[2];
	res[3] = _x[3] + pp[3];
	res[4] = _x[4] + pp[4];
	res[5] = _x[5] + pp[5];
	return( res );
    }

    ParticlePCyl operator-( const ParticlePCyl &pp ) const { 
	ParticlePCyl res;
	res[0] = _x[0] - pp[0];
	res[1] = _x[1] - pp[1];
	res[2] = _x[2] - pp[2];
	res[3] = _x[3] - pp[3];
	res[4] = _x[4] - pp[4];
	res[5] = _x[5] - pp[5];
	return( res );
    }

    ParticlePCyl operator*( double x ) const { 
	ParticlePCyl res;
	res[0] = _x[0]*x;
	res[1] = _x[1]*x;
	res[2] = _x[2]*x;
	res[3] = _x[3]*x;
	res[4] = _x[4]*x;
	res[5] = _x[5]*x;
	return( res );
    }

    friend ParticlePCyl operator*( double x, const ParticlePCyl &pp );
};


inline std::ostream &operator<<( std::ostream &os, const ParticlePCyl &pp )
{
    os << std::setw(12) << pp(0) << " "
       << std::setw(12) << pp(1) << " "
       << std::setw(12) << pp(2) << " "
       << std::setw(12) << pp(3) << " "
       << std::setw(12) << pp(4) << " "
       << std::setw(12) << pp(5);
    return( os );
}


inline ParticlePCyl operator*( double x, const ParticlePCyl &pp )
{
    ParticlePCyl res;
    res[0] = pp[0]*x;
    res[1] = pp[1]*x;
    res[2] = pp[2]*x;
    res[3] = pp[3]*x;
    res[4] = pp[4]*x;
    res[5] = pp[5]*x;
    return( res );
}


/*! \brief %Particle point class for 3D.
 *
 *  A collection of particle coordinates (t,x,vx,y,vy,z,vz)to describe
 *  time, location and velocity in 3D.
 */
class ParticleP3D : public ParticlePBase 
{
    double         _x[7];   /*!< \brief Point coordinates */

public:

    /*! \brief Default constuctor.
     */
    ParticleP3D() {}

    /*! \brief Constructor for 3D particle point.
     */
    ParticleP3D( double t, double x, double vx, double y, double vy, double z, double vz ) {
	_x[0] = t; _x[1] = x; _x[2] = vx; _x[3] = y; _x[4] = vy; _x[5] = z; _x[6] = vz;
    }

    /*! \brief Returns geometry mode.
     */
    static geom_mode_e geom_mode() { return(MODE_3D); }

    /*! \brief Returns number of dimensions for geometry.
     */
    static size_t dim() { return(3); }

    /*! \brief Returns number of coordinates used for particle point.
     */
    static size_t size() { return(7); }

    /*! \brief Returns time derivatives \a dxdt of coordinates at time
     *  \a t and coordinates \a x = (x,vx,y,vy,z,vz) for one particle.
     *
     *  The calculation of particle trajectory is done by integrating
     *  the Lorentz equation in a form of a set of ordinary
     *  differential equations. In the case of 3D the set is:
     *  \f{eqnarray*} \frac{dx}{dt}   &=& v_x \\ 
                      \frac{dy}{dt}   &=& v_y \\ 
                      \frac{dz}{dt}   &=& v_z \\ 
                      \frac{dv_x}{dt} &=& a_x = \frac{q}{m}(E_x + v_y B_z - v_z B_y) \\ 
                      \frac{dv_y}{dt} &=& a_y = \frac{q}{m}(E_y + v_z B_x - v_x B_z) \\
                      \frac{dv_z}{dt} &=& a_z = \frac{q}{m}(E_z + v_x B_y - v_y B_x) \f}
     */
    static int get_derivatives( double t, const double *x, double *dxdt, void *data );

    /*! \brief Return the number of trajectory intersections with
     *  plane \a crd = val on the trajectory from \a x1 to \a
     *  x2. Intersection points are appended to vector \a intsc.
     */
    static int trajectory_intersections_at_plane( std::vector<ParticleP3D> &intsc, 
						  int crd, double val,
						  const ParticleP3D &x1, const ParticleP3D &x2 );

    /*! \brief Returns the location of particle point in Vec3D.
     */
    Vec3D location() const { return( Vec3D( _x[1], _x[3], _x[5] ) ); }

    /*! \brief Returns the velocity of particle point in Vec3D.
     */
    Vec3D velocity() const { return( Vec3D( _x[2], _x[4], _x[6] ) ); }

    /*! \brief Returns speed of particle.
     */
    double speed() { return( sqrt(_x[2]*_x[2] + _x[4]*_x[4] + _x[6]*_x[6]) ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator[]( int i ) { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator[]( int i ) const { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator()( int i ) { return( _x[i] ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator()( int i ) const { return( _x[i] ); }

    ParticleP3D operator+( const ParticleP3D &pp ) const { 
	ParticleP3D res;
	res[0] = _x[0] + pp[0];
	res[1] = _x[1] + pp[1];
	res[2] = _x[2] + pp[2];
	res[3] = _x[3] + pp[3];
	res[4] = _x[4] + pp[4];
	res[5] = _x[5] + pp[5];
	res[6] = _x[6] + pp[6];
	return( res );
    }

    ParticleP3D operator-( const ParticleP3D &pp ) const { 
	ParticleP3D res;
	res[0] = _x[0] - pp[0];
	res[1] = _x[1] - pp[1];
	res[2] = _x[2] - pp[2];
	res[3] = _x[3] - pp[3];
	res[4] = _x[4] - pp[4];
	res[5] = _x[5] - pp[5];
	res[6] = _x[6] - pp[6];
	return( res );
    }

    ParticleP3D operator*( double x ) const { 
	ParticleP3D res;
	res[0] = _x[0]*x;
	res[1] = _x[1]*x;
	res[2] = _x[2]*x;
	res[3] = _x[3]*x;
	res[4] = _x[4]*x;
	res[5] = _x[5]*x;
	res[6] = _x[6]*x;
	return( res );
    }

    friend ParticleP3D operator*( double x, const ParticleP3D &pp );
};


inline std::ostream &operator<<( std::ostream &os, const ParticleP3D &pp )
{
    os << std::setw(12) << pp(0) << " "
       << std::setw(12) << pp(1) << " "
       << std::setw(12) << pp(2) << " "
       << std::setw(12) << pp(3) << " "
       << std::setw(12) << pp(4) << " "
       << std::setw(12) << pp(5) << " "
       << std::setw(12) << pp(6);
    return( os );
}


inline ParticleP3D operator*( double x, const ParticleP3D &pp )
{
    ParticleP3D res;
    res[0] = pp[0]*x;
    res[1] = pp[1]*x;
    res[2] = pp[2]*x;
    res[3] = pp[3]*x;
    res[4] = pp[4]*x;
    res[5] = pp[5]*x;
    res[6] = pp[6]*x;
    return( res );
}




/* ************************************************************************************************ *
 * Particle classes                                                                                 *
 * ************************************************************************************************ */


/*! \brief %Particle base class
 *
 *  A parent class for common (non-templated) parts of %Particle classes.
 */
class ParticleBase
{
protected:

    particle_status_e   _status;       /*!< \brief Status of particle */
    double              _IQ;           /*!< \brief Current or charge of particle
					*
					*   In time-dependent simulations particles are localized and IQ 
					*   is charge, in time-independent simulation IQ is current.
					*   In case of 2D simulations this is in units A/m or C/m.
					*   In 3D and cylindrically symmetric simulations unit is A or C.
					*/
    double              _qm;           /*!< \brief Ratio q/m [C/kg] */

    ParticleBase( double IQ, double q, double m ) 
	: _status(PARTICLE_OK), _IQ(IQ), _qm(q/m) {}

    ~ParticleBase() {}

public:

    /*! \brief Return particle status.
     */
    particle_status_e get_status() { return( _status ); }

    /*! \brief Set particle status.
     */
    void set_status( particle_status_e status ) { _status = status; }

    /*! \brief Return current or charge of particle [A/C].
     */
    double IQ() const { return( _IQ ); }

    /*! \brief Return q/m ratio [C/kg].
     */
    double qm() const { return( _qm ); }

};


/*! \brief %Particle class in some geometry.
 *
 *  Contains a templated description of one particle in some
 *  geometry. The geometry is selected throught templation, where PP
 *  is either ParticleP2D, ParticlePCyl or ParticleP3D. The particle
 *  classes for defined geometry modes can be used through defined
 *  types Particle2D, ParticleCyl and Particle3D.
 */
template<class PP> class Particle : public ParticleBase
{
    std::vector<PP>     _trajectory;   /*!< \brief Container of trajectory points */
    PP                  _x;            /*!< \brief Current coordinates of particle */

public:

    /*! \brief Constructor for particle.
     *
     *  Make new particle with 
     *  \param IQ Current or charge of the macroparticle (A/C).
     *  \param q  charge state of microscopic particle (C).
     *  \param m  mass of microscopic particle (kg).
     *  \param x  particle coordinates (m and m/s).
     */
    Particle( double IQ, double q, double m, const PP &x ) 
	: ParticleBase(IQ,q,m), _x(x) {}

    /*! \brief Destructor.
     */
    ~Particle() {}

    /*! \brief Operator for pointing to coordinate data.
     */
    double &operator()( int i ) { return( _x(i) ); }

    /*! \brief Operator for pointing to coordinate data.
     */
    const double &operator()( int i ) const { return( _x(i) ); }

    /*! \brief Returns the location of particle in Vec3D.
     */
    Vec3D location() const { return( _x.location() ); }

    /*! \brief Returns the velocity of particle in Vec3D.
     */
    Vec3D velocity() const { return( _x.velocity() ); }

    /*! \brief Return reference to coordinate data.
     */
    PP &x() { return( _x ); }

    /*! \brief Return const reference to coordinate data.
     */
    const PP &x() const { return( _x ); }

    /*! \brief Return reference to trajectory data.
     */
    PP &traj( int i ) { return( _trajectory[i] ); }

    /*! \brief Return const reference to trajectory data.
     */
    const PP &traj( int i ) const { return( _trajectory[i] ); }

    /*! \brief Return number of trajectory points of particle.
     */
    size_t traj_size( void ) const { return( _trajectory.size() ); }

    /*! \brief Add trajectory point to the end of the trajectory.
     */
    void add_trajectory_point( const PP &x ) { _trajectory.push_back( x ); }

    /*! \brief Define trajectory by copying.
     */
    void copy_trajectory( const std::vector<PP> &traj ) { _trajectory = traj; }

    /*! \brief Clears the particle trajectory.
     */
    void clear_trajectory( void ) { _trajectory.clear(); }

    /*! \brief Prints internal data to std::cout.
     */
    void debug_print( void ) const {
	size_t a;
	std::cout << "**Particle\n";
	if( _status == PARTICLE_OK )
	    std::cout << "status = PARTICLE_OK\n";
	else if( _status == PARTICLE_OUT )
	    std::cout << "status = PARTICLE_OUT\n";
	else if( _status == PARTICLE_COLL )
	    std::cout << "status = PARTICLE_COLL\n";
	else if( _status == PARTICLE_BADDEF )
	    std::cout << "status = PARTICLE_BADDEF\n";
	else if( _status == PARTICLE_TIME )
	    std::cout << "status = PARTICLE_TIME\n";
	else if( _status == PARTICLE_NSTP )
	    std::cout << "status = PARTICLE_NSTP\n";
	else
	    std::cout << "status = Unknown\n";
	std::cout << "IQ = " << _IQ << "\n";
	std::cout << "q/m = " << _qm << "\n";
	std::cout << "x = (";
	for( a = 0; a < _x.size()-1; a++ )
	    std::cout << _x(a) << ", ";
	std::cout << _x(a) << ")\n";
	std::cout << "Trajectory:\n";
	for( a = 0; a < _trajectory.size(); a++ ) {
	    std::cout << "x[" << a << "] = (";
	    uint32_t b;
	    for( b = 0; b < _trajectory[a].size()-1; b++ )
		std::cout << _trajectory[a](b) << ", ";
	    std::cout << _trajectory[a](b) << ")\n";
	}
    }
};


/*! \brief Particle class in 2D.
 *
 *  A typedef for templated class Particle.
 */
typedef Particle<ParticleP2D>  Particle2D;


/*! \brief Particle class in Cylindrical symmetry.
 *
 *  A typedef for templated class Particle.
 */
typedef Particle<ParticlePCyl> ParticleCyl;


/*! \brief Particle class in 3D.
 *
 *  A typedef for templated class Particle.
 */
typedef Particle<ParticleP3D>  Particle3D;



/*! \brief Temporary data bundle for particle iterators.
 */
struct ParticleIteratorData {
    ScalarField          *_scharge;       /*!< \brief Space charge field or NULL. */
    const Efield         *_efield;        /*!< \brief Electric field or NULL. */
    const VectorField    *_bfield;        /*!< \brief Magnetic field or NULL. */
    const Geometry       *_g;             /*!< \brief Geometry. */
    double                _qm;            /*!< \brief Precalculated q/m. */

    ParticleIteratorData( ScalarField *scharge, const Efield *efield, 
			  const VectorField *bfield, const Geometry *g )
	: _scharge(scharge), _efield(efield), _bfield(bfield), 
	  _g(g), _qm(0.0) {}
};



#endif













