/* *****************************************************************************
** *****************************************************************************
**
** CurveNDimension Lib is a C++ library for handling n-dimensional curves, or,
** actually, n-dimensional polygonal lines.
**
** Copyright (C) 2011-2015  Emilio Vital Brazil - emilio.brazil@gmail.com and
**                          Claudio Esperanca   - esperanc@cos.ufrj.br
**
** This file is part of CurveNDimension.
**
** CurveNDimension is free software: you can redistribute it and/or modify
** it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as published by
** the Free Software Foundation, either version 2.1 of the License, or
** (at your option) any later version.
**
** CurveNDimension is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with CurveNDimension.  If not, see <http://www.gnu.org/licenses/>.
**
** *****************************************************************************
** ****************************************************************************/

#ifndef VECTORN
#define VECTORN 1


// Define the tolerance to compare two vectors
#ifndef ZERO_VECTOR
    #define ZERO_VECTOR 1e-15
#endif

#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>

/// A N-dimensional Vector class
template < typename real = double , unsigned int dim = 3 >
class VectorN {
protected:
    real coord[dim]; //!< Where the coordinates are actually stored.

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
    
    /// \brief Coordinate indexing
    /// \param i  coordinate index.
    /// \return  reference to the i'th coordinate value.
    inline real& operator[] (unsigned i) { assert(i<dim); return coord[i]; }

    /// \brief Coordinate indexing
    /// \param i  coordinate index.
    /// \return  value of the i'th coordinate value.
    inline real  operator[] (unsigned i) const { assert(i<dim); return coord[i]; }

    /// \brief Coordinate indexing
    /// \param i  coordinate index.
    /// \return  value of the i'th coordinate value.
    inline real  at(unsigned i) const { assert(i<dim); return coord[i]; }

    //
    // Access by coordinate names
    //

    /// \brief Returns a reference to the first element.
    inline real& x(void) {
        return this->coord[0];
    }

    /// \brief Returns the value of the first coordinate.
    inline real x(void) const {
        return this->coord[0];
    }

    /// \brief The y coordinate is usually the second coordinate.
    ///
    /// Returns a reference to the second coordinate. Must have
    /// at least dimension 2.
    inline real& y(void) {
        assert (dim>1);
        return this->coord[1];
    }

    /// \brief Returns a reference to the second coordinate.
    ///
    /// Must have at least dimension 2.
    inline real y(void) const {
        assert (dim>1);
        return this->coord[1];
    }


    /// \brief The z coordinate is usually the third coordinate.
    ///
    /// Returns a reference to the third coordinate. Must have
    /// at least dimension 3.
    inline real& z(void) {
        assert (dim>2);
        return this->coord[2];
    }

    /// \brief Returns the value of the third coordinate.
    ///
    /// Must have at least dimension 3.
    inline real z(void) const {
        assert (dim>2);
        return this->coord[2];
    }

    //
    // Binary arithmetic methods
    //
    
    /// \brief Vector sum
    /// \param v another vector
    /// \return sum vector
    inline VectorN<real,dim> operator+(const VectorN<real,dim>& v) const {
        VectorN<real,dim> r;
        for (unsigned i = 0; i < dim; i++) r.coord[i] = coord[i]+v.coord[i];
        return r;
    }

    /// \brief Vector sum
    /// \param alpha a real value
    /// \return sum vector
    inline VectorN<real,dim> operator+(real alpha) const {
        VectorN<real,dim> r;
        for (unsigned i = 0; i < dim; i++) r.coord[i] = coord[i]+alpha;
        return r;
    }

    /// \brief Vector difference
    /// \param v another vector
    /// \return difference vector
    inline VectorN<real,dim> operator-(const VectorN<real,dim>& v) const {
        VectorN<real,dim> r;
        for (unsigned i = 0; i < dim; i++) r.coord[i] = coord[i]-v.coord[i];
        return r;
    }

    /// \brief Vector difference
    /// \param alpha a real value
    /// \return difference vector
    inline VectorN<real,dim> operator-(real alpha) const {
        VectorN<real,dim> r;
        for (unsigned i = 0; i < dim; i++) r.coord[i] = coord[i]-alpha;
        return r;
    }

    /// \brief Multiplication by scalar
    /// \param s scalar
    /// \return product vector
    template<typename scalar>
    inline VectorN<real,dim> operator*(scalar s) const {
        VectorN<real,dim> r;
        for (unsigned i = 0; i < dim; i++) r.coord[i] = coord[i] * (real) s;
        return r;
    }
    
    /// \brief Multiplication by inverse scalar
    /// \param s scalar
    /// \return  product vector
    inline VectorN<real,dim> operator/(real s) const { 
        VectorN<real,dim> r;
        for (unsigned i = 0; i < dim; i++) r.coord[i] = coord[i] / s;
        return r;
    }


    /// \brief Symmetric operator (unary minus)
    /// \return : symmetric vector
    inline VectorN<real,dim> operator-() const {
        VectorN<real,dim> r = *this;
        for (unsigned i = 0; i < dim; i++) r.coord[i] = -r.coord[i];
        return r;
    } 
    
    /// \brief Dot product
    /// \param v: another vector
    /// \return : dot product
    inline real operator*(const VectorN<real,dim>& v) const {
        real sum = 0;
        for (unsigned i = 0; i < dim; i++) sum += coord[i]*v.coord[i];
        return sum;
    }

    /// \brief Dot product
    /// \param v: another vector
    /// \return : dot product
    inline real dot(const VectorN<real,dim>& v) const {
        real sum = 0;
        for (unsigned i = 0; i < dim; i++) sum += coord[i]*v.coord[i];
        return sum;
    }
    
    /// Returns the square of the norm of the vector
    inline real norm2( void ) const { return (*this)*(*this); }

    /// Returns the Euclidian norm of the vector
    inline real norm( void ) const { return (real) sqrt((double) this->norm2()); }

    //
    // Assignment and in-place arithmetic methods
    //

    /// \brief Assignment
    /// \param v another vector
    /// \return reference to this (altered) vector
    inline VectorN<real,dim>& operator=(const VectorN<real,dim>& v) {
        for (unsigned i = 0; i < dim; i++) coord[i] = v.coord[i];
        return *this;
    }

    /// \brief Vector sum assignment
    /// \param v another vector
    /// \return reference to this (altered) vector
    inline VectorN<real,dim>& operator+=(const VectorN<real,dim>& v) {
        for (unsigned i = 0; i < dim; i++) coord[i] += v.coord[i];
        return *this;
    }

    /// \brief Vector sum assignment
    /// \param alpha a real value
    /// \return reference to this (altered) vector,
    ///         where alpha was summed in all coordinates
    inline VectorN<real,dim>& operator+=(real alpha) {
        for (unsigned i = 0; i < dim; i++) coord[i] += alpha;
        return *this;
    }

    /// \brief Vector difference assignment
    /// \param v another vector
    /// \return reference to this (altered) vector
    inline VectorN<real,dim>& operator-=(const VectorN<real,dim>& v) {
        for (unsigned i = 0; i < dim; i++) coord[i] -= v.coord[i];
        return *this;
    }

    /// \brief Vector difference assignment
    /// \param alpha a real value
    /// \return reference to this (altered) vector
    inline VectorN<real,dim>& operator-=(real alpha) {
        for (unsigned i = 0; i < dim; i++) coord[i] -= alpha;
        return *this;
    }


    /// \brief Multiplication by scalar assignment
    /// \param s scalar
    /// \return  reference to this (altered) vector
    inline VectorN<real,dim>& operator*=(real s) {
        for (unsigned i = 0; i < dim; i++) coord[i] *= s;
        return *this;
    }

    /// \brief Multiplication by inverse scalar assignment
    /// \param s scalar
    /// \return reference to this (altered) vector
    inline VectorN<real,dim>& operator/=(real s) { return (*this) *= (1/s); }

    //
    // Comparison operators
    //

    /// \brief Equality operator
    /// \param v another vector
    /// \return whether this and v have equal coordinates within the error margin
    bool operator==(const VectorN<real,dim> &v) const {
        for (unsigned i = 0; i < dim; i++) {
            real d = std::abs(this->coord[i]-v.coord[i]);
            if (d > ZERO_VECTOR ) return false;
        }
        return true;
    }

    /// \brief Not Equal operator
    /// \param v another vector
    /// \return whether this and v have at least one distinct coordinate within
    ///   the error margin
    bool operator!=(const VectorN<real,dim>& v) const {
        return !(*this == v);
    }

    /// \brief less than operator using lexicographical order
    /// \param v another vector
    bool operator<(const VectorN<real,dim>& v) const {
        unsigned i = 0 ;
        for (i = 0; i < dim; i++)
            if (std::abs(this->coord[i]-v.coord[i]) > ZERO_VECTOR ) break ;

        i = ( i < dim ) ? i : i-1 ;
        return ( this->coord[i] < v.coord[i] );
    }

    /// \brief less or equal than operator
    /// \param v another vector
    /// \return true if all entries of vector are less
    /// than or equal to the
    ///         respctive entry of v
    bool operator<=(const VectorN<real,dim>& v) const {
        for (unsigned i = 0; i < dim; i++) {
            if (this->coord[i] > v.coord[i]) return false;
        }
        return true;
    }

    /// \brief less than operator
    /// \param alpha a real value
    /// \return true if all entries of vector are less than the alpha
    bool operator<(real alpha) const {
        for (unsigned i = 0; i < dim; i++) {
            if (this->coord[i] >= alpha ) return false;
        }
        return true;
    }

    /// \brief less than operator
    /// \param alpha a real value
    /// \return true if all entries of vector are less or then than the alpha
    bool operator<=(real alpha) const {
        for (unsigned i = 0; i < dim; i++) {
            if (this->coord[i] > alpha ) return false;
        }
        return true;
    }

    /// \brief greater than operator using lexicographical order
    /// \param v another vector
    bool operator>(const VectorN<real,dim>& v) const {
        unsigned i = 0 ;
        for (i = 0; i < dim; i++)
            if (std::abs(this->coord[i]-v.coord[i]) > ZERO_VECTOR ) break;
        i = ( i < dim ) ? i : i - 1 ;
        return ( this->coord[i] <= v.coord[i] );
    }

    /// \brief greater than or equal to operator
    /// \param v another vector
    /// \return true if all entries of vector are greater than or equal to the
    ///         respctive entry of v
    bool operator>=(const VectorN<real,dim>& v) const {
        for (unsigned i = 0; i < dim; i++) {
            if (this->coord[i] < v.coord[i]) return false;
        }
        return true;
    }

    /// \brief greater than operator
    /// \param alpha a real value
    /// \return true if all entries of vector are greater than to the alpha
    bool operator>(real alpha) const {
        for (unsigned i = 0; i < dim; i++) {
            if (this->coord[i] <= alpha ) return false;
        }
        return true;
    }

    /// \brief less than or equal to operator
    /// \param alpha a real value
    /// \return true if all entries of vector are greater than or equal to the alpha
    bool operator>=(real alpha) const {
        for (unsigned i = 0; i < dim; i++) {
            if (this->coord[i] < alpha ) return false;
        }
        return true;
    }


};


/// Makes scalar multiplication commutative.
template < typename real, unsigned int dim >
inline VectorN<real,dim> operator*(real s, const VectorN<real,dim>& v) {
    return v*s;
}

/// Makes scalar vector sum commutative.
template < typename real, unsigned int dim >
inline VectorN<real,dim> operator+(real s, const VectorN<real,dim>& v) {
    return v+s;
}

/// Makes scalar vector difference commutative.
template < typename real, unsigned int dim >
inline VectorN<real,dim> operator-(real s, const VectorN<real,dim>& v) {
    return -(v-s);
}

//
// Primitive function definitions
//

template < typename real >
inline real minReal(real a , real b)
{
    return ( a < b ) ? a : b ;
}

template < typename real >
inline real maxReal(real a , real b)
{
    return ( a > b ) ? a : b ;
}

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

/// \brief Prints a VectorN onto an output stream
/// \param out: an output stream
/// \param v : a VectorN
/// \return: the modified output stream
template < typename real, unsigned int dim >
inline std::ostream& operator<<(std::ostream& out, const VectorN<real,dim>& v)
{
    for (unsigned i = 0; i < dim; i++) {
        out << v [i] << " ";
    }
    return out;
}

/// \brief Reads a VectorN from an input stream
/// \param in: an output stream
/// \param v : a VectorN
/// \return: the modified input stream
template < typename real, unsigned int dim >
inline std::istream& operator>>(std::istream& in, VectorN<real,dim>& v)
{
    for (unsigned i = 0; i < dim; i++) {
        in >> v [i];
    }
    return in;
}


/// A N-dimensional Point class
template < typename real = double , unsigned int dim = 3 >
class PointN : public VectorN < real, dim > {
public:
    //
    // Standard constructors
    //
    
    /// Empty constructor
    PointN () : VectorN<real,dim> () {}
    
    /// Constructor from an array of coordinates
    PointN (const real v[]) : VectorN<real,dim> (v) {}
    
    /// Constructor from another VectorN
    PointN (const VectorN<real,dim>& v) : VectorN<real,dim> (v) {}

    /// Constructor from a single coordinate 
    PointN ( real v ) : VectorN<real,dim> () {
        for (unsigned int i = 0; i < dim; i++) this->coord [i] = v;
    }

    /// Constructor from two coordinates - must be bidimensional
    PointN ( real x , real y) : VectorN<real,dim> () {
        assert (dim == 2);
        this->coord [0] = x;
        this->coord [1] = y;
    }

    /// Constructor from three coordinates - must be tridimensional
    PointN ( real x,  real y, real z) : VectorN<real,dim> () {
        assert (dim == 3);
        this->coord [0] = x;
        this->coord [1] = y;
        this->coord [2] = z;
    }

    //
    // Auxiliaary Functions
    //

    /// \brief Set coordinates from a vector
    /// \param v: a vector
    inline void set( const VectorN<real,dim>& v ) { 
        for ( unsigned int i = 0; i < dim; i++) this->coord [i] = v[i];
    }
    
    /// Casting to Vector
    inline VectorN<real,dim> toVector ( void ) { return VectorN<real,dim> (*this); }
    
//    /// Affine multiplication by a scalar
//    inline PointN<real,dim> affineSum( real alpha , const PointN<real,dim>& p) const {
//        return (PointN<real,dim>) (alpha * p);
//    }
    
    /// \brief Get the minimum coordinates.
    /// \param p another point.
    /// \return point whose coordinates are the minimum between this point's coordinates
    ///    and p's
    inline PointN<real,dim> min(const PointN<real,dim> & p) const {
        PointN<real,dim> r;
        for ( unsigned int i = 0; i < dim; i++) r[i] = minReal(this->at(i),p[i]);
        return r;
    }

    /// \brief Get maximum coordinates
    /// \param p another point.
    /// \return point whose coordinates are the maximum between this point's coordinates
    ///     and p's
    inline PointN<real,dim> max(const PointN<real,dim> & p) const {
        PointN<real,dim> r;
        for ( unsigned int i = 0; i < dim; i++) r[i] = maxReal(this->at(i),p[i]);
        return r;
    }
    
    /// Squared distance between two points
    inline real dist2 (const PointN<real,dim> & p) const {
        PointN<real,dim> q(*this);
        return (p-q).norm2();
    }
    
    /// Euclidian distance between two points
    inline real dist (const PointN<real,dim> & p) const {
        return ((*this) - p).norm();
    }

    
};

#endif

   
