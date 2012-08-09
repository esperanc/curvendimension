#ifndef __CURVEN_HPP__
#define __CURVEN_HPP__

#include <vector>
#include <queue>
#include <cassert>
#include <limits>
#include "vectorn.hpp"
#include "catmullrom.hpp"

/// Infinity constant
#ifndef INF
#define INF std::numeric_limits<real>::max()
#endif
    
/// Epsilon constant
#ifndef EPS
#define EPS std::numeric_limits<real>::epsilon()
#endif

/**
* @brief N-dimensional polygonal chain.
*
* This templated class represents an open or closed Polygonal Chain in n dimensions. 
* It requires Point and Vector types for which a standard implementation is 
* provided. If desired, the user may provide a custom implementation of these 
* classes or, better yet, derive and customize them.
* 
* @param real is the type used for storing each cartesian coordinate of each point.
* @param dim is the number of dimensions of each point.
* @param Point each vertex of the polygonal line is stored as an object of this class.
* @param Vector objects of this class are returned whenever a vector in dim dimensions is
*       required.
*
**/
template < typename real , unsigned int dim , class Point = PointN<real,dim> , 
           class Vector = VectorN<real,dim> >
class PolygonalCurve {

protected:
	/// @brief a Smooth Step function (see smoothstep in wikipedia).
	///
	/// A simple smooth function that maps the unit interval to itself. 
	/// Smoothness is guaranteed by having the derivatives at 0 and 1 to
	/// be zero. (see smoothstep in Wikipedia)
	///
	/// @param x a value that is constrained to be between 0 and 1 
	/// @return  a value between 0 and 1
	static inline double smoothStep (double x) {
		if (x < 0) x = 0;
		else if (x>1) x = 1;
	    return x*x*(3 - 2*x);
	    // alternatively
	    // return x*x*x*(x*(x*6 - 15) + 10);
	}

    /// @brief class used by the douglasPeuckerRank method.
    ///
    /// This is used as an element in the priority queue needed
    /// by the douglasPeuckerRank method.
    struct DPitem {
        unsigned first;    ///< index of the first vertex of the range
        unsigned last;     ///< index of the last vertex of the range
        unsigned farthest;  ///< index between first and last which is farthest from
                            /// the line segment defined by the vertices first and last
        double dist;        ///< squared perpendicular distance from the line segment
        
        /// @brief Constructor.
        ///
        /// Builds an item corresponding to the a vertex range of a polygonal curve.
        ///
        /// @param f  index of first point in poly
        /// @param l  index of last point in poly
        /// @param poly polygonal curve which is being generalized
        DPitem (unsigned f, unsigned l, const PolygonalCurve<real,dim,Point,Vector>& poly) 
            : first(f), last(l) 
        {
            assert (last-first > 1); 
            // Compute farthest point and its distance to line first-last
            const Point& p0 = poly [first];
            const Point& p1 = poly [last]; 
            Vector v = p1-p0; // Direction vector from p0 to p1
            unitize (v);      // Make it a unit vector
            farthest = first;
            dist = -1;
            for (unsigned i = first+1; i < last; i++) {
                const Point& p = poly [i];
                Point pr = p0+((p-p0)*v)*v;  // p projected onto line p0-p1
                double d = pr.dist2(p); // Squared distance
                if (d > dist) {  // Keep index of the farthest point
                    farthest = i;
                    dist = d;
                }
            }   
        }
        
        /// Operator <
        inline bool operator < (const DPitem& other) const { return dist < other.dist; }
    }; 

public:    
    /// Empty constructor - builds an open polygonal line with 0 vertices.
    PolygonalCurve( void ) { pIsClosed = false ; }

    /// @brief Copy constructor.
    /// 
    /// Builds the curve as a copy of another polygonal curve.
    /// @param curve: another polygonal curve.
    PolygonalCurve( const PolygonalCurve<real,dim,Point,Vector> &curve ) { 
        this->copy(curve) ; 
    }

    /// @brief Constructor from an array.
    ///
    /// Constructor from an array (stl vector data type) of Points.
    /// @param points an array of points.
    PolygonalCurve( const std::vector<Point> &points ) {
        pPoints = points;
        if ( pPoints.size() != 0 && points[ pPoints.size()-1 ] == pPoints[ 0 ] )
        {
            pIsClosed = true ;
            pPoints.pop_back();
        }
        else pIsClosed = false ;
    }

    /// @brief Copy helper method.
    ///
    /// Replaces this object with a copy of another curve. 
    /// @param curve Makes this object a copy of curve.
    virtual void copy( const PolygonalCurve<real,dim,Point,Vector> &curve ) {
        pPoints = curve.pPoints;
        pIsClosed = curve.pIsClosed;
    }

    /// @brief Assignment operator
    ///
    /// Makes this object a copy of the given curve. 
    /// @param curve: another polygonal curve.
    /// @return
    virtual PolygonalCurve &operator=( const PolygonalCurve<real,dim,Point,Vector> &curve ) {
        this->copy(curve);
        return (*this);
    }

    /// @brief Indexing operator.
    /// 
    /// Returns the reference to point at index i
    /// @param i the index.
    virtual Point& operator[](unsigned i)       { assert(i<size()) ; return pPoints[i] ; }
    
    /// @brief Indexing operator.
    /// 
    /// Returns the value of point at index i
    /// @param i: the index.
    virtual Point  operator[](unsigned i) const { assert(i<size()) ; return pPoints[i] ; }

    /// Number of points of the curve
    virtual unsigned int size( void ) const { return pPoints.size() ; }


    /// @brief Replaces the points of this curve by points in an array.
    ///
    /// Changes the internal representation of the geometry to points in a given array
    /// @param points: an array of points
    virtual void setPoints( const std::vector<Point> &points ) { pPoints = points ; }
    
    /// @brief Returns a copy of the geometry as an array of points
    ///
    /// @return an array of points.
    virtual std::vector<Point> getPoints( void ) const { return pPoints ; }

    /// @brief Resets line to an empty open curve.
    virtual void clear( void ) { pPoints.clear() ; pIsClosed = false; }

    /// @brief Returns a copy of the i'th point in the curve.
    /// @param i: index.
    virtual Point   at( unsigned int i ) const {
        assert(i<size()) ;
        return pPoints[i] ;
    }
    
    /// @brief Returns a copy of the first point in the curve. 
    /// 
    /// Requires a non-empty curve.
    /// @return the first point.
    virtual Point   atBegin( void ) const { assert(size()>0) ; return pPoints[0] ; }
    
    /// @brief Returns the last point in the curve. 
    /// 
    /// Requires a non-empty curve.
    /// @return the last point.
    virtual Point   atEnd( void ) const {  assert(size()>0) ; return pPoints[ size()-1 ] ; }

    /// @brief Computes C(u) for u in [0,1].
    ///
    /// Computes C(u), assuming the curve is parameterized by arc length,
    /// i.e., C(0) is the first point, C(1) is the last point.
    /// @param u a value between 0 and 1;
    /// @param p (output) the point at C(u);
    /// @param alpha (output) if C(u) is between vertex i and vertex i+1, returns
    ///     the interpolation ratio between them, such that alpha=0 means C(u) is P(i)
    ///     and alpha=1 means C(u) = P(i+1)
    /// @return the index of the point just before u is reached. 
    virtual unsigned int  eval( real u , Point& p , real& alpha ) const {
        if( size() == 0 ) return 0;
        real l = length();
        if( u >= 1.0 )
        {
            p = at(size()-1);
            return size()-1;
        }
        if( u <= 0.0 || l == 0 )
        {
            p = at(0);
            alpha = 0.0;
            return 0;
        }
        real il = 1.0/l;
        real aux1 = 0 ;
        real aux2 = 0 ;
        real lenTmp = 0 ;
        unsigned int index = 0 ;
        Vector vtmp;
        while( aux1 <= u )
        {
            lenTmp = aux1;
            index++;
            vtmp = pPoints[index]-pPoints[index-1];
            aux2 = vtmp.norm()*il;
            aux1 += aux2;
        }
        alpha = (u - lenTmp)/aux2;
        index--;
        p = at(index) + alpha*vtmp;
        return index;
    }

    /// @brief Computes C(u) for u in [0,1].
    ///
    /// Computes C(u), assuming the curve is parameterized by arc length,
    /// i.e., C(0) is the first point, C(1) is the last point.
    /// @param u a value between 0 and 1;
    /// @return C(u). 
    virtual Point eval( real u ) const { Point p ; eval( u , p ) ; return p ; }
    
    /// @brief Computes C(u) for u in [0,1].
    ///
    /// Computes C(u), assuming the curve is parameterized by arc length,
    /// i.e., C(0) is the first point, C(1) is the last point.
    /// @param u a value between 0 and 1;
    /// @param p (output) the point at C(u);
    /// @return the index of the point just before u is reached. 
    virtual unsigned int  eval( real u , Point& p ) const { real alpha ; return eval( u , p , alpha ); }
    
    /// @brief Element assignment.
    /// Alters the i'th point
    /// @param index the index.
    /// @param p  new point.
    virtual void setPoint( unsigned int index , const Point& p ){ 
        assert(index<size()) ; pPoints[index] = p ; 
    }
    
    /// @brief Curve extension.
    ///
    /// Adds another point to the curve.
    /// @param p  new point.
    virtual void add( const Point& p ) { 
        pPoints.push_back(p) ;
    }

    /// @brief Point insertion.
    ///
    /// Inserts a point at a given index.
    /// @param index the index.
    /// @param p point to be inserted.
    virtual void insert( unsigned int index , const Point& p ) { 
        if( index < size()) pPoints.insert( pPoints.begin()+index , p ); else add(p) ; 
    }


    /// @brief Curve trimming.
    ///
    /// Removes the last point of the curve
    virtual void pop_back(void) { pPoints.pop_back(); }
    
    /// @brief Curve extension.
    ///
    /// Appends another point to the end of the curve.
    /// @param p the point to be appended.
    virtual void push_back( const Point& p ) { pPoints.push_back(p) ; }

    /// @brief Closes or opens curve.
    ///
    /// Makes this a closed/open curve.
    /// @param b whether the curve must be closed (true) or not (false).
    virtual void close( bool b = true ) { pIsClosed = b ; }
    
    /// Makes this an open curve (polyline).
    virtual void open( void ) { pIsClosed = false ; }
    
    /// @brief Tells whether the curve is closed.
    /// @return true iff curve is closed.
    virtual bool isClosed() const { return pIsClosed ; }

    /// @brief Chaikin simplification of the polygonal line.
    /// @param numberOfSimplifications number of simplifications
    virtual void chaikinFilter( unsigned int numberOfSimplifications = 1 ) {
        if (size() == 0) return;
        if( isClosed() )
        {
            for( unsigned int n = 0 ; n < numberOfSimplifications ; ++n )
            {
                std::vector<Point> tmp;
                // make sure we do not lose any points!!
                if ( size() <= 6 ) return;

                // step over every 2 points
                for( unsigned int i = 1 ; i < ( size() - 2 ) ; i += 2 )
                {
                    // get original points
                    const Vector p0 = pPoints[i-1].toVector();
                    const Vector p1 = pPoints[i  ].toVector();
                    const Vector p2 = pPoints[i+1].toVector();
                    const Vector p3 = pPoints[i+2].toVector();

                    // calculate the original point
                    Vector Q = -0.25f*p0 + 0.75f*p1 + 0.75f*p2 - 0.25f*p3;

                    // add to new curve
                    tmp.push_back( Point(Q) );
                }
                if( size()%2 != 0 )
                {
                    unsigned int i = size() - 1 ;
                    // get original points
                    const Vector p0 = pPoints[i-1].toVector();
                    const Vector p1 = pPoints[i  ].toVector();
                    const Vector p2 = pPoints[0].toVector();
                    const Vector p3 = pPoints[1].toVector();
                    // calculate the original point
                    Vector Q = -0.25f*p0 + 0.75f*p1 + 0.75f*p2 - 0.25f*p3;

                    // add to new curve
                    tmp.push_back( Point(Q) );
                }
                else
                {
                    unsigned int i = size() - 2 ;
                    // get original points
                    const Vector p0 = pPoints[i-1].toVector();
                    const Vector p1 = pPoints[i  ].toVector();
                    const Vector p2 = pPoints[i+1].toVector();
                    const Vector p3 = pPoints[0].toVector();
                    // calculate the original point
                    Vector Q = -0.25f*p0 + 0.75f*p1 + 0.75f*p2 - 0.25f*p3;

                    // add to new curve
                    tmp.push_back( Point(Q) );

                }
                // copy over pPoints
                pPoints = tmp ;
            }
        }
        else
        {
            for( unsigned int n = 0 ; n < numberOfSimplifications ; ++n )
            {
                std::vector<Point>  tmp;
                // make sure we do not loose any points!!
                if (size() <= 8) return;

                // keep the first point
                tmp.push_back( pPoints[0] );
                Point p2 = -0.5f*pPoints[0].toVector() + pPoints[1].toVector() +
                        0.75f*pPoints[2].toVector() - 0.25f*pPoints[3].toVector();
                tmp.push_back( Point(p2) );

                // step over every 2 points
                for( unsigned int i = 2 ; i < ( size() - 5 ) ; i += 2 )
                {
                    // get original points
                    const Vector p0 = pPoints[i    ].toVector();
                    const Vector p1 = pPoints[i + 1].toVector();
                    const Vector p2 = pPoints[i + 2].toVector();
                    const Vector p3 = pPoints[i + 3].toVector();
                    // calculate the original point
                    Vector Q = -0.25f*p0 + 0.75f*p1 + 0.75f*p2 - 0.25f*p3 ;
                    // add to new curve
                    tmp.push_back( Point(Q) );
                }
                unsigned int lastIndex = size() - 1 ;
                Vector pL = -0.25f*pPoints[lastIndex-3].toVector() + 
                        0.75f*pPoints[lastIndex-2].toVector() +
                        pPoints[lastIndex-1].toVector() - 0.50f*pPoints[lastIndex].toVector();
                tmp.push_back( Point(pL) );
                tmp.push_back(pPoints[ lastIndex ]);
                // copy over pPoints
                pPoints = tmp ;
            }
        }
    }

    /// @brief Chaikin supersampling.
    /// @param numberOfsimplifications number of subdivisions.
    virtual void chaikinSubDivide( unsigned int numberOfsimplifications = 1 ){
        if (size() == 0) return;
        if( isClosed() )
        {
            for( unsigned int n = 0 ; n < numberOfsimplifications ; ++n )
            {
                std::vector<Point>  tmp;
                // step over every 2 points
                for( unsigned int i = 1 ; i < size() ; ++i )
                {
                    // get original points
                    const Vector p0 = pPoints[i-1].toVector();
                    const Vector p1 = pPoints[i  ].toVector();

                    // calculate the original point
                    Vector Q = 0.75f*p0 + 0.25f*p1 ;

                    tmp.push_back( Point(Q) );

                    Vector R ;
                    R = 0.25f*p0 + 0.75f*p1 ;
                    tmp.push_back( Point(R) );
                }
                // get original points
                const Vector p0 = atEnd().toVector();
                const Vector p1 = atBegin().toVector();

                // calculate the original point
                Vector Q = 0.75f*p0 + 0.25f*p1 ;

                tmp.push_back( Point(Q) );

                Vector R ;
                R = 0.25f*p0 + 0.75f*p1 ;
                tmp.push_back( Point(R) );
                // copy over pPoints
                pPoints = tmp ;
            }
        }
        else
        {
            for( unsigned int n = 0 ; n < numberOfsimplifications ; ++n )
            {
                std::vector<Point>  tmp;
                // keep the first point
                tmp.push_back( pPoints[0] );

                // step over every 2 points
                for( unsigned int i = 1 ; i < size() ; ++i )
                {
                    // get original points
                    const Vector p0 = pPoints[i-1].toVector();
                    const Vector p1 = pPoints[i  ].toVector();

                    // calculate the original point
                    Vector Q = 0.75f*p0 + 0.25f*p1 ;
                    tmp.push_back( Point(Q) );
                    Vector R ;
                    R = 0.25f*p0 + 0.75f*p1 ;
                    tmp.push_back( Point(R) );
                }
                tmp.push_back(pPoints[ size()-1 ] );
                // copy over pPoints
                pPoints = tmp ;
            }
        }
    }

    /// @brief Simple supersampling.
    /// @param step makes sure each two consecutive vertices are separated by at least this
    ///  amount.
    virtual void superSample( real step ){
        if (size() == 0) return;

        PolygonalCurve tmpLine;
        tmpLine.add( pPoints[0] );
        for( unsigned int i = 0 ; i <  size() - 1 ; ++i )
        {
            Vector tmpV = pPoints[i+1] -pPoints[i] ;
            real tmpNorm = tmpV.norm() ;
            if( tmpNorm > step )
            {
                Vector base = pPoints[i].toVector() ;
                real invNorm = 1/tmpNorm ;
                Vector unitV = invNorm*tmpV ;
                unsigned int j = 1 ;
                while( (j*step) <= ( tmpNorm - step ) )
                {
                    Vector newP = base + (j*step)*unitV ;
                    tmpLine.add( Point( newP ) );
                    ++j;
                }
            }
            tmpV = pPoints[i+1] -tmpLine.atEnd() ;
            tmpNorm = tmpV.norm() ;
            if( tmpNorm > 0.8*step )
            {
                tmpLine.add( pPoints[i+1] );
            }
        }
        if(isClosed())
        {
            Vector tmpV = atBegin() - atEnd() ;
            real tmpNorm = norm( tmpV ) ;
            //        real step = 0.5f;
            if( tmpNorm > step )
            {
                Vector base = atEnd().toVector() ;
                real invNorm = 1/tmpNorm ;
                Vector unitV = invNorm*tmpV ;
                unsigned int j = 1 ;
                while( (j*step) <= ( tmpNorm - step ) )
                {
                    Vector newP = base + (j*step)*unitV ;
                    tmpLine.add( Point( newP ) );
                    ++j;
                }
            }
        }
        tmpLine.close( isClosed() );
        *this = tmpLine;
    }

    /// @brief Applies the superSample(step) and chaikinFilter( numberOfInteractions ).
    ///
    /// Useful to clean up hand-drawn input curves.
    /// @param step: maximum distance between consecutive vertices in supersampling.
    /// @param numberOfInteractions: how many simplification steps are performed.
    virtual void lineFilter( real step = 0.5 , unsigned int numberOfInteractions = 5 ){
        if (size() == 0) return;
        this->superSample( step );
        this->chaikinFilter( numberOfInteractions );
    }
    
    /// Mean filter, i.e, each vertex v(i) at i, is replaced by (v(i-1)+v(i)*3+v(i+1))/5.
    virtual void meanFilter( void ){
        if( size() > 1 )
        {
            Vector v;
            std::vector<Point> tmp;
            if(isClosed())
            {
                v = ( atEnd().toVector()+3*at(0).toVector()+at(1).toVector() )/5.0;
                tmp.push_back( Point(v[0],v[1],v[2]) );
            }
            else tmp.push_back( atBegin() );

            for( unsigned int i = 1 ; i < size()-1 ; ++i )
            {
                v = ( at(i-1).toVector()+3*at(i).toVector()+at(i+1).toVector() )/5.0;
                tmp.push_back( Point(v[0],v[1],v[2]) );
            }
            if (isClosed())
            {
                v = ( at(size()-2).toVector()+3*atEnd().toVector()+at(0).toVector() )/5.0;
                tmp.push_back( Point(v[0],v[1],v[2]) );
            }
            else tmp.push_back( atEnd() );
            pPoints = tmp;
        }
    }

    /// @brief edits curve by dragging a point to a new position.
    ///
    /// Deforms the curve by displacing vertex p[k] to a new position pk. 
	/// Neighboring vertices whose distance along the curve to p[k] are less than
	/// | pk-p[k] | * factor are also displaced smoothly. 
	/// @param k index of the anchor point which is the basis for the deformation.
	/// @param pk New position of point k.
	/// @param factor proportionality constant for determining the neighborhood of 
	///   the anchor point - larger values mean a bigger neighborhood.
	inline void smoothDeform (unsigned k, const Point& pk, double factor = 1.0) {
		Vector v = pk - at(k);
		double vlen = v.norm();
		double maxLen = vlen * factor;
		double d = 0;
		Point p0 = (*this) [k];
		(*this) [k] = pk;
		Point prev = p0;
		for (int i = k-1; i >= 0; i--) {
			Point pi = at(i);
			double dincr = (pi-at(i+1)).norm();
			d += (pi-prev).norm();
			if (d > maxLen) break;
			prev = pi;
			pi += v * smoothStep((maxLen-d)/maxLen);
			(*this) [i] = pi;
		}
		d = 0;
		prev = p0;
		for (int i = k+1; i < size(); i++) {
			Point pi = at(i);
			d += (pi-prev).norm();
			if (d > maxLen) break;
			prev = pi;
			pi += v * smoothStep((maxLen-d)/maxLen);
			(*this) [i] = pi;
		}
	}

    /// @brief Catmull-Rom spline interpolation.
    ///
    /// Creates a Catmull-Rom interpolated curve with
    /// the vertices of this Polygonal Curve as control points.
    /// @param result (output): where the interpolated curve is to be stored - any previous
    ///     content is erased.
    /// @param maxdist distance between interpolated points should be not greater than this.
    inline void catmull (PolygonalCurve<real,dim,Point,Vector>& result, double maxdist = 1.0) {
    
        result.clear(); // Make sure result has no vertices
        double d2 = maxdist*maxdist; // Use quicker squared distance 
        
		// Catmullrom interpolator
		CatmullRomBlend<real, dim> crb;

        // How many segments?
        int nsegments = isClosed()? size() : size()-1;
        
        // Emit first point
        result.add (at(0));
        
        // Emit at least nsegments points
        int i = 0;
        while (i < nsegments) {
        
            // The two middle control points
            Point p1 = at(i % size());
            Point p2 = at((i+1) % size());
                      
            if (p1.dist2(p2) > d2) {
                // Perform subdivision
                
                // The two extremity control points
                Point p0 = i == 0 ? (isClosed() ? atEnd() : at(0)) : at (i-1);
                Point p3 = i+2 < size() || isClosed() ? at((i+2) % size()) : atEnd();

				// The anchor point against which to measure distance
				Point anchor = p1;
				
				// The u value corresponding to the anchor
				double uanchor = 0.0;
				
				// A stack of u's (parameter values)
				std::vector<real> ustack;
				ustack.push_back (1.0);
				ustack.push_back (0.5);
				
				// Perform iterative subdivision
				while (ustack.size()>0) {
					real u = ustack.back();
					ustack.pop_back();
					Point p;
                    crb.blendPoint (u,&p0.x(),&p1.x(),&p2.x(),&p3.x(),&p.x());
					if (p.dist2(anchor) <= d2) {
						result.add (p);
						anchor = p;
						uanchor = u;
					} else {
						real unew = (u+uanchor) / 2;
						ustack.push_back (u);
						ustack.push_back (unew);
					}
				}
            } 
            else {
				// Add segment point
            	result.add (p2);
			}
            
            i++;
        }
    }
    
    
    /// @brief Performs a Douglas-Peucker analysis of the curve.
    ///
    /// Returns an array of ranks of vertices according to the
    /// generalization order imposed by the Douglas Peucker algorithm.
    /// Thus, if the i'th element has value k, then vertex i would be the (k+1)'th
    /// element to be included in a generalization (simplification) of this polyline.
    /// @param tol do not consider vertices farther than this value. Disconsidered
    ///     vertices are marked with -1 in the result.
    /// @return  a rank array.
    std::vector<int> douglasPeuckerRank (double tol) const {
        
        // A priority queue of intervals to subdivide
        std::priority_queue<DPitem> pq;
        
        // The result vector
        std::vector<int> r (size(), -1);
        
        // Put the first and last vertices in the result
        r [0] = 0;
        r [size()-1] = 1;
        
        // Add first interval to pq
        if (size()>1) pq.push (DPitem (0,size()-1, (*this)));
        
        // The rank counter 
        int rank = 2;
        
        // Recursively subdivide up to tol
        while (!pq.empty()) {
            DPitem item = pq.top();
            pq.pop();
            if (item.dist < tol*tol) break; // All remaining points are closer 
            r [item.farthest] = rank++;
            if (item.farthest > item.first+1) 
                pq.push (DPitem (item.first, item.farthest, (*this)));
            if (item.last > item.farthest+1) 
                pq.push (DPitem(item.farthest, item.last, (*this))); 
        }
        
        return r;
    }                         
    
    /// @brief Douglas-Peucker generalization.
    ///   
    /// Returns a simplification (generalization) of this polygonal line
    /// using the Douglas Peucker algorithm with the given tolerance.
    /// @param result where the result curve is stored (previous contents are erased).
    /// @param tol remove vertices farther than this value.
    void douglasPeuckerSimplify (PolygonalCurve<real,dim,Point,Vector>& result, double tol) const {
        assert (size()>0);
        std::vector<int> rank = douglasPeuckerRank(tol);
        result.clear();
        for (unsigned i = 0; i < size(); i++) {
            if (rank[i] >= 0) { 
                result.add(at(i));
            }
        }
    }
    
    /// @brief Douglas-Peucker generalization.
    ///
    /// Returns a simplification (generalization) of this polygonal line
    /// with at most n vertices using the Douglas Peucker algorithm.
    /// @param result where the result curve is stored (previous contents are erased).
    /// @param n maximum number of points in the result.
    void douglasPeuckerDecimate (PolygonalCurve<real,dim,Point,Vector>& result, unsigned n) const {
        assert (size()>0);
        assert (n >= 2 && n <= size());
        std::vector<int> rank = douglasPeuckerRank(0.0);
        result.clear();
        for (unsigned i = 0; i < size(); i++) {
            if (rank[i] < n) { 
                result.add(at(i));
            }
        }
    }

    /// Makes line follow the reverse circulation.
    virtual void reverse( void ){
        std::vector<Point> tmp;
        for( unsigned int i = size() ; i > 0 ; --i )
        {
            tmp.push_back( pPoints[i-1] );
        }
        pPoints = tmp;
    }

    /// @brief Joins another curve with this one.
    ///
    /// Joins this curve with an additional curve, testing if the first and last points
    /// are not repeated, i.e., they are nearer than err.
    /// @param l: another curve.
    /// @param err: minimum distance to consider two points as the same.
    virtual void join( const PolygonalCurve& l , real err = EPS ){
        if(l.size() == 0 ) return;
        real dist = INF ;
        if( this->size() > 0 ) dist = ( atEnd() - l.atBegin()).norm();
        unsigned int start = ( dist < err ) ? 1 : 0 ;
        for( unsigned int i = start ; i < l.size()-1 ; ++i ) this->add(l.at(i));
        dist = (l.atEnd() - atBegin()).norm() ;
        if( dist < err ) close( true );
        else add(l.atEnd());
    }

    /// @brief Tangent at u.
    ///
    /// Estimates a tangent vector at a given point C(u), where u in [0,1].
    /// @param u: relative distance along the curve from the beginning of the curve.
    /// @return : the estimated tangent at C(u).
    virtual Vector tangentEval( real u ) const {
        Vector tan;
        if( size() > 1 )
        {
            Point p;
            real alpha;
            unsigned int index = eval( u , p , alpha );
            if( alpha > EPS )
            {
                if( index != size()-1 ) tan = at(index+1) - at(index);
                else if ( isClosed() ) tan = at(0) - at(index);
                else tan = at(index-1) - at(index);
            }
            else tan = tangent(index);
        }
        return tan;
    }

    /// @brief Estimated tangent at a given vertex. 
    /// @param index: index between 0 and size()-1.
    /// @return estimated tangent.
    virtual Vector tangent( unsigned int index ) const {
        Vector tan;
        if( size() > 1 )
        {
            if( index < size()-1 && index > 0 ) tan = ( at(index+1) - at(index-1) )*.5 ;
            else if( isClosed() ) tan = ( at( (index+1)%size() ) - at( (int(index) - 1) % size() ) ) *.5;
            else if ( index > 0 ) tan = at(index) - at(index-1);
            else tan = at(1) - at(0);
        }
        return tan;
    }

    /// @brief Tangent at u.
    ///
    /// Smoothly interpolated tangent at a point along the curve.
    /// @param u: relative distance along the curve from the beginning of the curve.
    /// @return : estimated tangent.
    virtual Vector tangentEvalContinuous( real u ) const {
        Vector tan;
        if( size() > 1 )
        {
            Point p;
            real alpha;
            unsigned int index = eval( u , p , alpha );
            if( index != size()-1 ) tan = alpha * tangent(index) + (1.0-alpha)* tangent(index+1);
            else  tan = tangent(index);
        }
        return tan;
    }

    /// @brief Bounding box of the curve.
    /// @param min (output) point with minimum coordinates.
    /// @param max (output) point with maximum coordinates.
    virtual void minMax( Point& min , Point& max ) const {
        max = -INF, min = INF ;

        for( unsigned int i = 0 ; i < size() ; i++ )
        {
            for( unsigned int j = 0 ; j < dim ; j++ )
            {
                max[j] = ( at(i)[j] > max[j] ) ? at(i)[j] : max[j] ;
                min[j] = ( at(i)[j] < min[j] ) ? at(i)[j] : min[j] ;
            }
        }

    }

    /// Computes the arc length of the curve.
    /// @return the length of the curve.
    virtual real length( void ) const {
        real len = 0;
        Vector vtmp;
        for( unsigned int i = 1 ; i < size() ; i++ )
        {
            vtmp = pPoints[i]-pPoints[i-1];
            len += vtmp.norm();
        }
        if(isClosed()) vtmp = atBegin() - atEnd() , len += vtmp.norm() ;
        return len;
    }
    
    /// @brief Computes the length of the curve up to a given vertex.
    /// @param k index of the vertex (between 0 and size()-1)
    /// @return length.
    virtual real length( unsigned int k ) const {
        real len = 0;

        k = ( k < size() )? k+1 : size();

        for( unsigned int i = 1 ; i < k ; i++ )
        {
            Vector vtmp = pPoints[i]-pPoints[i-1];
            len += vtmp.norm();
        }
        return len;
    }

    /// @brief Computes the relative distance between a vertex and the beginning of the curve.
    /// @param i index of the vertex (between 0 and size()-1)
    /// @return  relative length (between 0 and 1)
    virtual real getParameter( unsigned int i ) const { return length(i)/length() ; }

    /// @brief Adds a vector to all points in the curve.
    /// @param v displacement vector.
    virtual void translate ( const Vector& v ) { 
        for( unsigned int i = 0 ; i < size() ; ++i ) pPoints[i] += v;
    }

    /// @brief Computes the distance between a point and the curve.
    /// @param p a point.
    /// @return  Euclidean distance between p and curve.
    virtual real distanceTo( const Point& p ) const { 
        unsigned int index ; real u ; Point pr ; 
        return projectPoint( p , index , u , pr ) ; 
    }
    
    /// @brief Computes distance between a point and the curve.
    /// @param p input point.
    /// @param index (output) index of the vertex just before the segment on which the 
    ///     projection of p on the curve lies.
    /// @param parameterU (output) relative length of the projected point with respect to the total
    ///     length of the curve.
    /// @param pr (output) closest point on the curve to p.
    /// @return  distance between p and pr.
    virtual real projectPoint( const Point& p , unsigned int& index , real& parameterU , Point& pr ) const {
        real dist = INF, uTmp1 = 0.0, uTmp2 = 0.0 ;
        parameterU = 0.0;
        if( size() > 0 )
        {
            index = 0 ;
            Point prTmp = at(0);
            pr = prTmp;
            dist = ( atBegin() - p ).norm() ;
            real distTMP = INF , invL = 0.0 ;
            if( length()>EPS) invL = 1.0/length();
            for( unsigned int i = 1 ; i < size() ; i++ )
            {
                uTmp1 = uTmp2;
                unsigned int iTmp = i;
                Vector u = at(i) - at(i-1);
                Vector v = p - at(i-1);
                Vector w = p - at(i);
                real normQ = u.norm2();
                real alpha = 0.0;
                real distU = u.norm();
                real distV = v.norm();
                real distW = w.norm();

                if( normQ > EPS ) alpha = (u*v)/normQ ;

                if( alpha > 0.0 && alpha < 1.0 )
                {
                    prTmp = at(i-1) + alpha*u;
                    distTMP = (prTmp-p).norm();
                    uTmp1 += distU*invL*alpha;
                }
                else if ( distV < distW ) distTMP =  distV, prTmp = at(i-1) ;
                else distTMP = distW , ++iTmp , uTmp1 += distU*invL , prTmp = at(i) ;

                if( distTMP < dist ) dist = distTMP , index = iTmp - 1 , parameterU = uTmp1 , pr = prTmp;
                uTmp2 += distU*invL;
            }
            if(isClosed())
            {
                uTmp1 = uTmp2;
                unsigned int iTmp =  size();
                Vector u = atBegin() - atEnd();
                Vector v = p - atEnd();
                Vector w = p - atBegin();
                real normQ = u.norm2();
                real alpha = 0.0;
                real distU = u.norm();
                real distV = v.norm();
                real distW = w.norm();

                if( normQ > EPS ) alpha = (u*v)/normQ ;

                if( alpha > 0.0 && alpha < 1.0 )
                {
                    prTmp = atEnd() + alpha*u;
                    distTMP = (prTmp-p).norm();
                    uTmp1 += distU*invL*alpha;
                }
                else if ( distV < distW ) distTMP =  distV, prTmp = atEnd() ;
                else distTMP = distW , iTmp = 0 , uTmp1 += distU*invL , prTmp = atBegin() ;

                if( distTMP < dist ) dist = distTMP , index = iTmp - 1 , parameterU = uTmp1 , pr = prTmp;
            }
        }
        return dist;
    }

    /// @brief Performs a scale operation on the curve.
    /// @param v the diagonal of the scale transformation matrix.
    virtual void scale(const Vector & v) {
        for( unsigned int i = 0 ; i < size() ; i++ )
            for( unsigned int j = 0 ; j < dim ; j++ ) pPoints[i][j]*=v[j];
    }

    /// Returns the centroid (barycenter) of the vertex points
    virtual Point centroid( void ) const {
        Vector barTmp;

        for( unsigned int i = 0 ; i < size() ; i++ ) barTmp += at(i).toVector();

        barTmp /= size() ;

        return Point( barTmp );
    }


    /// @brief Splits this curve at the given index. 
    /// 
    /// This curve will be trimmed to the part 
    /// before index and the method returns the portion of the original curve after index.
    /// The point at(index) is placed in both.
    /// @param index index of split point.
    /// @return cut out curve segment.
    virtual PolygonalCurve<real,dim,Point,Vector> split( unsigned int index ) {

        PolygonalCurve<real,dim,Point,Vector> newCurve;

        if(index<size())
        {
            std::vector<Point> tmp ;
            tmp = std::vector<Point>( pPoints.begin()+index , pPoints.end() ) ;
            newCurve.setPoints(tmp) ;
            tmp = std::vector<Point>( pPoints.begin(),pPoints.begin()+index+1 ) ;
            setPoints(tmp);
            if(isClosed()) newCurve.add( atBegin() );
            open();
            newCurve.open();
        }

        return newCurve;
    }
    
    /// @brief Split curve at two vertices.
    ///
    /// Creates 2 or 3 new curves by splitting this curve. If this curve is closed,
    /// 2 curves (c1 and c2) are returned, and, if not, 3 (c1, c2 and c3).
    /// The new curves are this curve split at the indices index1 and index2.
    /// @param index0 index of first split point.
    /// @param index1 index of second split point.
    /// @param c1 (output) First segment.
    /// @param c2 (output) Second segment.
    /// @param c3 (output) Third segment.    
    virtual bool split( unsigned int index0 , unsigned int index1 ,
        PolygonalCurve<real,dim,Point,Vector>& c1 ,
        PolygonalCurve<real,dim,Point,Vector>& c2 , 
        PolygonalCurve<real,dim,Point,Vector>& c3  ) const 
    {
        bool b = index1<size() && index0<size() ;
        if( b )
        {
            std::vector<Point> tmp ;
            if(isClosed())
            {
                tmp = std::vector<Point>( pPoints.begin() , pPoints.begin()+index0+1 ) ;
                c1.setPoints(tmp) ;
                tmp = std::vector<Point>( pPoints.begin()+index0 , pPoints.begin()+index1+1 ) ;
                c2.setPoints(tmp) ;
                tmp = std::vector<Point>( pPoints.begin()+index1 , pPoints.end() ) ;
                c3.setPoints(tmp) ;
                c3.join(c1);
                c1 = c3;
                c3.clear();
            } else
            {
                tmp = std::vector<Point>( pPoints.begin() , pPoints.begin()+index0+1 ) ;
                c1.setPoints(tmp) ;
                tmp = std::vector<Point>( pPoints.begin()+index0 , pPoints.begin()+index1+1 ) ;
                c2.setPoints(tmp) ;
                tmp = std::vector<Point>( pPoints.begin()+index1 , pPoints.end() ) ;
                c3.setPoints(tmp) ;
            }
        }
        return b;
    }
    
    /// @brief Splits curve at a vertex.
    /// 
    /// Creates 2 new curves (c1 and c2) by splitting the current curve.
    /// The new curves are this curve splitted at the given index. The
    /// original curve is not modified.
    ///
    /// @param index index where the curve will be split.
    /// @param c1 (output) first segment.
    /// @param c2 (output) second segment.
    /// @return true iff c2 is not empty.
    virtual bool split( unsigned int index , PolygonalCurve<real,dim,Point,Vector>& c1 ,
                        PolygonalCurve<real,dim,Point,Vector>& c2 ) const {
        c1 = (*this);
        c2 = c1.split( index );
        return c2.size() > 0;
    }
    
    /// @brief splits curve at the closest point to a given point.
    ///
    /// Creates 2 new curves (c1 and c2) by splitting the current curve.
    /// The new curves are this curve splitted at the projection of p.
    /// If addPoint is true, p is added to c1 and c2 instead of its projection.
    ///
    /// @param p the split point is the curve point closest to this.
    /// @param c1 (output) first segment.
    /// @param c2 (output) second segment.
    /// @param addPoint if true, add p rather than its projection to both segments.
    /// @param err do not add any point if it is closer than this distance to a curve vertex.
    /// @return true iff c2 is not empty.
    virtual bool split( const Point& p , 
        PolygonalCurve<real,dim,Point,Vector>& c1 ,
        PolygonalCurve<real,dim,Point,Vector>& c2 , 
        bool addPoint = false , 
        real err = EPS ) const 
    {
        real u;
        unsigned int index;
        Point pr;

        c1 = (*this);

        c1.projectPoint( p , index , u , pr );
        index += (index>0) ? 0 : 1 ;
        if(addPoint)
        {
            if( (c1.at(index)-p).norm() > err && (c1.at(index-1)-p).norm() > err )
                ++index , c1.insert(index,p);
        }
        else if( (c1.at(index)-pr).norm() > err && (c1.at(index-1)-pr).norm() > err )
        {
            ++index ;
            std::cerr << " index = " << index << std::endl;
            c1.insert( index , pr );
        }
        c2 = c1.split(index) ;
        return c2.size() > 0;
    }

    /// @brief Edits a curve segment with another.
    ///
    /// Returns a new curve that is the result of using this curve as an over-sketch curve.
    /// The difference between the returned curve and the original curve is 
    /// saved in rest; this curve may have its orientation changed.
    /// The parameter "err" controls if two points are the same, and aproximationFactor 
    /// is how many times the start and or 
    /// points should be far from the curve to be interpreted as an extension 
    /// rather than an over-sketch.
    /// @param curve (input) the curve to be edited.
    /// @param rest (out) altered part of the original curve.
    /// @param err two points are considered the same if closer than this distance.
    /// @param aproximationFactor determines how close a point must be to be considered an oversketch.
    /// @result the edited curve.
    virtual PolygonalCurve<real,dim,Point,Vector> overSketch( const 
           PolygonalCurve<real,dim,Point,Vector>& curve ,
           PolygonalCurve<real,dim,Point,Vector>& rest ,
           real err = EPS , 
           real aproximationFactor = 64.0 ) {
        Point p0 = atBegin();
        Point p1 = atEnd();

        Point pr0 , pr1 , prTmp;
        real u0, u1, uTmp, distTmp;
        unsigned int index0 , index1, indexTmp;

        real dist0 = curve.projectPoint( p0 , index0 , u0 , pr0 ) ;
        real dist1 = curve.projectPoint( p1 , index1 , u1 , pr1 ) ;

        if( u0 > u1 )
        {
            prTmp = pr0;
            uTmp = u0;
            indexTmp = index0;
            distTmp = dist0;
            pr0 = pr1;
            u0 = u1;
            index0 = index1;
            dist0 = dist1;
            pr1 = prTmp;
            u1 = uTmp;
            index1 = indexTmp;
            dist1 = distTmp;
            this->reverse();
        }

        PolygonalCurve<real,dim,Point,Vector> A , B , C , curveTmp0 , curveTmp1  ;

        if( dist0 > aproximationFactor*err && dist1 < aproximationFactor*err )
        {
            curve.split( pr1 , rest , curveTmp0 , false , err );
            this->join( curveTmp0 , err );
            curveTmp0 = (*this);
        }
        else if( dist0 < aproximationFactor*err && dist1 > aproximationFactor*err   )
        {
            curve.split( pr0 , curveTmp0 , rest , false , err );
            curveTmp0.join( (*this) , err ) ;
        }
        else
        {
            curve.split( pr0 , A , curveTmp0 , false , err );
            curveTmp0.split( pr1 , B , C , false , err );

            Point centOr = curve.centroid();

            curveTmp0 = A;
            curveTmp0.join((*this),err);
            curveTmp0.join(C,err);
            curveTmp0.close(isClosed());
            Point cent0 = curveTmp0.centroid();

            curveTmp1 = B;
            curveTmp1.reverse();
            curveTmp1.join((*this),err);
            curveTmp1.close();
            Point cent1 = curveTmp1.centroid();

            if( (centOr-cent0).norm2() < (centOr-cent1).norm2() ) rest = B ;
            else
            {
                curveTmp0 = curveTmp1;
                rest = C;
                rest.join(A,err);
            }
        }

        return curveTmp0 ;
    }

protected:
    std::vector<Point> pPoints;  ///< Where the points are actually stored.
    bool pIsClosed;              ///< Whether or not this is a closed curve.

};


#endif // __CURVEN_HPP__

