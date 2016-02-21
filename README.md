ntroduction

This is a C++ library for handling n-dimensional curves, or, actually, n-dimensional polygonal lines.

# Classes provided

* `template<typename real, unsigned int dim, class Point = PointN<real,dim>, class Vector = VectorN<real,dim>>`
    `class PolygonalCurve< real, dim, Point, Vector >`

    N-dimensional polygonal chain.

    This templated class represents an open or closed Polygonal Chain in n dimensions. It requires Point and Vector types for which a standard implementation is provided. If desired, the user may provide a custom implementation of these classes or, better yet, derive and customize them.

    Parameters:

    + `real`	is the type used for storing each cartesian coordinate of each point.

    + `dim`	is the number of dimensions of each point.

    + `Point`	each vertex of the polygonal line is stored as an object of this class.

    + `Vector`	objects of this class are returned whenever a vector in dim dimensions is required.
