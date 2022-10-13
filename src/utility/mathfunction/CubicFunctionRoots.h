#ifndef _CUBIC_FUNCTION_ROOTS_H_
#define _CUBIC_FUNCTION_ROOTS_H_

#include "../../main/MPM3D_MACRO.h"
//!> Calculate the roots of "a*x^3 + b*x^2 + c*x + d = 0" with Shengjin formulation
//!> Fan Shengjin. A new extracting formula and a new distinguishing means on the one variable cubic equation. pp. 91—98 .
inline void CubicFunctionRoots(MPM_FLOAT a, MPM_FLOAT b, MPM_FLOAT c, MPM_FLOAT d, Array3D& roots)
{
    MPM_FLOAT A = b*b - 3*a*c;
    MPM_FLOAT B = b*c - 9*a*d;
    MPM_FLOAT C = c*c - 3*b*d;
    MPM_FLOAT Delta = B*B - 4*A*C;

    if (fabs(Delta) <= MPM_EPSILON)
    {
        if (fabs(A) <= MPM_EPSILON || fabs(a) <= MPM_EPSILON)
        {
            roots[0] = roots[1] = roots[2] = 0;
            return;
        }

        MPM_FLOAT K = B/A;
        roots[0] = -b/a + K;
        roots[1] = roots[2] = -K/2;
    }
    else if (Delta < -MPM_EPSILON)
    {
        if (A < -MPM_EPSILON)
        {
            roots[0] = roots[1] = roots[2] = 0;
            return;
        }

        MPM_FLOAT T = (2*A*b - 3*a*B)/(2*pow(A, 1.5));
        MPM_FLOAT theta = acos(T);
        MPM_FLOAT ct = cos(theta/3);
        MPM_FLOAT st = sin(theta/3);
        MPM_FLOAT sqrtA = sqrt(A);
        roots[0] = (-b - 2*sqrtA*ct)/(3*a);
        roots[1] = (-b + sqrtA*(ct + SQRT3*st))/(3*a);
        roots[2] = (-b + sqrtA*(ct - SQRT3*st))/(3*a);
    }
    else
    {
        roots[0] = roots[1] = roots[2] = 0;
    }
}

#endif