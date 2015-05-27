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

#ifndef __CATMULLROM_HPP__
#define __CATMULLROM_HPP__
///
/// Cubic Catmull-Rom Interpolator.
///
/// @param real: scalar type to use (usually float or double)
/// @param dim: domain dimension
///
template <typename real, unsigned int dim = 3>
class CatmullRomBlend {
  
protected:

  real tau; ///< Tension
  
public:   
  
    /// @brief Constructor
    /// @param tension: the tension.
    CatmullRomBlend (real tension = 0.5) {
        tau = tension; 
    }

    /// @brief Given a parameter u, computes the blending factors for
    /// the four control points. 
    /// @param u: curve parameter. 
    /// @param factor (output): array of 4 blending factors.
    inline void blendFactors (real u, real factor []) const {
      real u1 = u;
      real u2 = u*u;
      real u3 = u2*u;
      factor [0] = -tau * u1 + 2 * tau * u2 - tau * u3;
      factor [1] = 1 + (tau-3) * u2 + (2 - tau) * u3;
      factor [2] = tau * u1 + (3 - 2*tau) * u2 + (tau - 2) * u3;
      factor [3] = -tau * u2 + tau * u3;
    }
  
    /// @brief Given a parameter u and the coordinates of 4 control points, returns 
    /// computes the interpolated point.
    /// @param u: curve parameter.
    /// @param p0: first control point (at least size dim).
    /// @param p1: second control point (at least size dim).
    /// @param p2: third control point (at least size dim).
    /// @param p3: fourth control point (at least size dim).
    /// @param p (output): interpolated point.
    inline void blendPoint (real u, const real p0[], const real p1[], 
                            const real p2[], const real p3[], real p[]) const {
        real bf [4];
        blendFactors (u,bf);
        for ( unsigned int j = 0 ; j <  dim ; j++) {
            p[j] = p0[j]*bf[0]+ 
                   p1[j]*bf[1]+
                   p2[j]*bf[2]+
                   p3[j]*bf[3];
        }
    }
};
#endif
