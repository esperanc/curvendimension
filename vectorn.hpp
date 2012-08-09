#ifndef VECTORN
#define VECTORN 1

#ifndef FEQ_EPS
#define FEQ_EPS 1e-6
#define FEQ_EPS2 1e-12
#endif

#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>

/// A N-dimensional Vector class
template < typename real = double , unsigned int dim = 3 >
class VectorN {
protected:
    real coord[dim]; ///< Where the coordinates are actually stored.

public:
    //
    // Standard constructors
    //
    
    /// Empty constructor
    VectorN () { 
        for (unsigned i = 0; i < dim; i++) coord[i] = 0; 
    }
    
    /// Constructor from an array of coordinates
    VectorN (const real v[]) {  
        for (unsigned i = 0; i < dim; i++) coord[i] = v[i]; 
    }
    
    /// Constructor from another VectorN
    VectorN (const VectorN<real,dim>& v) { 
        for (unsigned i = 0; i < dim; i++) coord[i] = v.coord[i]; 
    }
    
    //
    // Access methods
    //
    
    /// @brief Coordinate indexing 
    /// @param i  coordinate index.
    /// @return  reference to the i'th coordinate value.
    inline real& operator[] (unsigned i) { assert(i<dim); return coord[i]; }

    /// @brief Coordinate indexing 
    /// @param i  coordinate index.
    /// @return  value of the i'th coordinate value.
    inline real  operator[] (unsigned i) const { assert(i<dim); return coord[i]; }

    //
    // Comparison operators
    //
    
    /// @brief Equality operator
    /// @param v another vector
    /// @return whether this and v have equal coordinates within the error margin
    inline bool operator==(const VectorN<real,dim> &v) const {
        for (unsigned i = 0; i < dim; i++) {
            real d = coord[i]-v.coord[i]; 
            if (d > FEQ_EPS2) return false;
        }
        return true;
    }
    
    /// @brief Inequality operator
    /// @param v another vector
    /// @return whether this and v have at least one distinct coordinate within 
    ///   the error margin
    inline bool operator!=(const VectorN<real,dim>& v) const {
        return !(*this == v);
    }
    
    //
    // Assignment and in-place arithmetic methods
    //
    
    /// @brief Assignment
    /// @param v another vector
    /// @return reference to this (altered) vector
    inline VectorN<real,dim>& operator=(const VectorN<real,dim>& v) {
        for (unsigned i = 0; i < dim; i++) coord[i] = v.coord[i];
        return *this;
    };
    
    /// @brief Vector sum assignment
    /// @param v another vector
    /// @return reference to this (altered) vector
    inline VectorN<real,dim>& operator+=(const VectorN<real,dim>& v) {
        for (unsigned i = 0; i < dim; i++) coord[i] += v.coord[i];
        return *this;
    };
    
    /// @brief Vector difference assignment
    /// @param v another vector
    /// @return reference to this (altered) vector
    inline VectorN<real,dim>& operator-=(const VectorN<real,dim>& v) {
        for (unsigned i = 0; i < dim; i++) coord[i] -= v.coord[i];
        return *this;
    };
    
    /// @brief Multiplication by scalar assignment
    /// @param s scalar
    /// @return  reference to this (altered) vector
    inline VectorN<real,dim>& operator*=(real s) {
        for (unsigned i = 0; i < dim; i++) coord[i] *= s;
        return *this;
    };
    
    /// @brief Multiplication by inverse scalar assignment
    /// @param s scalar
    /// @return reference to this (altered) vector
    inline VectorN<real,dim>& operator/=(real s) { return (*this) *= (1/s); }

    //
    // Binary arithmetic methods
    //
    
    /// @brief Vector sum 
    /// @param v another vector
    /// @return sum vector
    inline VectorN<real,dim> operator+(const VectorN<real,dim>& v) const {
        VectorN<real,dim> r;
        for (unsigned i = 0; i < dim; i++) r.coord[i] = coord[i]+v.coord[i];
        return r;
    };
    
    /// @brief Vector difference
    /// @param v another vector
    /// @return difference vector
    inline VectorN<real,dim> operator-(const VectorN<real,dim>& v) const {
        VectorN<real,dim> r;
        for (unsigned i = 0; i < dim; i++) r.coord[i] = coord[i]-v.coord[i];
        return r;
    };
    
    /// @brief Multiplication by scalar
    /// @param s scalar
    /// @return product vector
    template<typename scalar>
    inline VectorN<real,dim> operator*(scalar s) const {
        VectorN<real,dim> r;
        for (unsigned i = 0; i < dim; i++) r.coord[i] = coord[i] * (real) s;
        return r;
    };
    
    /// @brief Multiplication by inverse scalar
    /// @param s scalar
    /// @return  product vector
    inline VectorN<real,dim> operator/(real s) const { 
        VectorN<real,dim> r;
        for (unsigned i = 0; i < dim; i++) r.coord[i] = coord[i] / s;
        return r;
    }


    /// @brief Symmetric operator (unary minus)
    /// @return : symmetric vector
    inline VectorN<real,dim> operator-() const {
        VectorN<real,dim> r = *this;
        for (unsigned i = 0; i < dim; i++) r.coord [i] = -r.coord[i];
        return r;
    } 
    
    /// @brief Dot product
    /// @param v: another vector
    /// @return : dot product
    inline real operator*(const VectorN<real,dim>& v) const {
        real sum = 0;
        for (unsigned i = 0; i < dim; i++) sum += coord[i]*v.coord[i];
        return sum;
    }
    
    /// Returns the square of the norm of the vector
    inline real norm2( void ) const { return (*this)*(*this); }

    /// Returns the Euclidian norm of the vector
    inline real norm( void ) const { return (real) sqrt((double) this->norm2()); }


    //
    // Access by coordinate names
    //
    
    /// @brief Returns a reference to the first element.
    inline real& x(void) {
        return this->coord[0];
    }

    /// @brief Returns the value of the first coordinate.
    inline real x(void) const {
        return this->coord[0];
    }

    /// @brief The y coordinate is usually the second coordinate.
    ///
    /// Returns a reference to the second coordinate. Must have
    /// at least dimension 2.
    inline real& y(void) {
        assert (dim>1);
        return this->coord[1];
    }

    /// @brief Returns a reference to the second coordinate. 
    ///
    /// Must have at least dimension 2.
    inline real y(void) const {
        assert (dim>1);
        return this->coord[1];
    }

    
    /// @brief The z coordinate is usually the third coordinate.
    ///
    /// Returns a reference to the third coordinate. Must have
    /// at least dimension 3.
    inline real& z(void) {
        assert (dim>2);
        return this->coord[2];
    }
    
    /// @brief Returns the value of the third coordinate. 
    ///
    /// Must have at least dimension 3.
    inline real z(void) const {
        assert (dim>2);
        return this->coord[2];
    }
};


/// Makes scalar multiplication commutative.
template < typename scalar, typename real, unsigned int dim >
inline VectorN<real,dim> operator*(scalar s, const VectorN<real,dim>& v) { 
    return v*(real)s; 
}

//
// Primitive function definitions
//

/// Returns the Euclidian norm of a Vector
template < typename real, unsigned int dim >
inline real norm(const VectorN<real,dim>& v)
{
    return v.norm();
}

/// Returns the squared length of a Vector
template < typename real, unsigned int dim >
inline real norm2(const VectorN<real,dim>& v)
{
    return v.norm2();
}

/// Makes a vector unit length
template < typename real, unsigned int dim >
inline real unitize(VectorN<real,dim>& v)
{
    double l=norm2(v);
    if( l!=1.0 && l!=0.0 )   {
    	l = sqrt(l);
    	v /= l;
    }
    return l;
}

/// @brief Prints a VectorN onto an output stream
/// @param out: an output stream
/// @param v : a VectorN
/// @return: the modified output stream
template < typename real, unsigned int dim >
inline std::ostream& operator<<(std::ostream& out, const VectorN<real,dim>& v)
{
    for (unsigned i = 0; i < dim; i++) {
        out << v [i] << " ";
    }
    return out;
}

/// @brief Reads a VectorN from an input stream
/// @param in: an output stream
/// @param v : a VectorN
/// @return: the modified input stream
template < typename real, unsigned int dim >
inline std::istream& operator>>(std::istream& in, VectorN<real,dim>& v)
{
    for (unsigned i = 0; i < dim; i++) {
        in >> v [i];
    }
    return in;
}


///
/// A N-dimensional Point class
///
template < typename real = double , unsigned int dim = 3 >
class PointN : public VectorN < real, dim > {
public:
    //
    // Standard constructors
    //
    
    /// Empty constructor
    PointN () : VectorN<real,dim> () {};
    
    /// Constructor from an array of coordinates
    PointN (const real v[]) : VectorN<real,dim> (v) {};
    
    /// Constructor from another VectorN
    PointN (const VectorN<real,dim>& v) : VectorN<real,dim> (v) {}; 

    /// Constructor from a single coordinate 
    PointN (const real v) : VectorN<real,dim> () {
        for (int i = 0; i < dim; i++) this->coord [i] = v;
    }

    /// Constructor from two coordinates - must be bidimensional
    PointN (const real x, const real y) : VectorN<real,dim> () {
        assert (dim == 2);
        this->coord [0] = x;
        this->coord [1] = y;
    }

    /// Constructor from three coordinates - must be tridimensional
    PointN (const real x, const real y, const real z) : VectorN<real,dim> () {
        assert (dim == 3);
        this->coord [0] = x;
        this->coord [1] = y;
        this->coord [2] = z;
    }

    /// @brief Set coordinates from a vector
    /// @param v: a vector
    inline void set( const VectorN<real,dim>& v ) { 
        for (int i = 0; i < dim; i++) this->coord [i] = v[i];
    }
    
    /// Casting to Vector
    inline VectorN<real,dim> toVector ( void ) { return VectorN<real,dim> (*this); }
    
    /// Affine multiplication by a scalar
    inline PointN<real,dim> affineSum( real alpha , const PointN<real,dim>& p) const {
        return (PointN<real,dim>) (alpha * p);
    }
    
    /// @brief Get point with minimum coordinates.
    /// @param p another point.
    /// @return point whose coordinates are the minimum between this point's coordinates
    ///    and p's
    inline PointN<real,dim> min(const PointN<real,dim> & p) const {
        PointN<real,dim> r;
        for (int i = 0; i < dim; i++) r[i] = min(*this[i],p[i]);
        return r;
    }

    /// @brief Get point with maximum coordinates 
    /// @param p another point.
    /// @return point whose coordinates are the maximum between this point's coordinates
    ///     and p's
    inline PointN<real,dim> max(const PointN<real,dim> & p) const {
        PointN<real,dim> r;
        for (int i = 0; i < dim; i++) r[i] = max(*this[i],p[i]);
        return r;
    }
    
    /// Squared distance between two points
    inline real dist2 (const PointN<real,dim> & p) const {
        return ((*this) - p).norm2();
    }
    
    /// Euclidian distance between two points
    inline real dist (const PointN<real,dim> & p) const {
        return ((*this) - p).norm();
    }

    
};

#endif

   
