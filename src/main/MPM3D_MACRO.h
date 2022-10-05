/*==============================================================
                            OpenMPM3D
    C-plus-plus code for 3-Dimensional Material Point Method
================================================================
    Copyright (C) 2022 - 

    Computational Dynamics Group
    Department of Engineering Mechanics
    School of Aerospace Engineering
    Tsinghua Univeristy
    Beijing 100084, P. R. China

    Corresponding Author: Xiong Zhang
    E-mail: xzhang@tsinghua.edu.cn
================================================================
    Info: Define common macros and parameters
    Code-writter: Ruichen Ni
    Date: 2022.9.28
==============================================================*/

#ifndef _MPM3D_MACRO_H_
#define _MPM3D_MACRO_H_

#include <string>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>
using namespace std;

#ifdef _MPM_DOUBLE
    typedef double MPM_FLOAT;
    const MPM_FLOAT MPM_EPSILON =   DBL_EPSILON;
    const MPM_FLOAT MPM_MIN =       DBL_MIN;
    const MPM_FLOAT MPM_MAX =       DBL_MAX;
    const MPM_FLOAT PI =            3.14159265358979323846;
    const MPM_FLOAT SQRT3 =         1.7320508075688773;
    constexpr MPM_FLOAT FourThird = 4.0/3.0;
#else
    typedef float MPM_FLOAT;
    const MPM_FLOAT MPM_EPSILON =   FLT_EPSILON;
    const MPM_FLOAT MPM_MIN =       FLT_MIN;
    const MPM_FLOAT MPM_MAX =       FLT_MAX;
    const MPM_FLOAT PI =            3.1415926536F;
    const MPM_FLOAT SQRT3 =         1.7320508076F;
    constexpr MPM_FLOAT FourThird = 4.0F/3.0F;
#endif

#ifdef _MPM_MASSIVE_PARTICLE
    typedef long MPM_STATS;
#else
    typedef int MPM_STATS;
#endif

typedef array<MPM_FLOAT, 3> Array3D;

namespace MPM{
    const string ProgramType = "MPM3D-CPP";
    const string MPM3DVersion = "1.";       //!< Release version
    const string FileVersion = "1.0";       //!< File format version
    const string MPM3DRevision = "0.0";     //!< Revision number
    const string CompileDate = __DATE__;
    const string CompileTime = __TIME__;

    enum MPMScheme
    {
        USL,        //!< update-stress-late formulation
        MUSL,       //!< modified USL formulation
        USF         //!< update-stress-first formulation
    };

    //!> Number of extra particle properties
    const MPM_STATS ExtraParticlePropertySum = 15;
    //!> Extra property list
    enum ExtraParticleProperty
    {
        SDxx, SDxy, SDxz,
              SDyy, SDyz,
                    SDzz,   //!< Deviatoric stress
        Exx, Exy, Exz,
             Eyy, Eyz,
                  Ezz,      //!< Strain
        epeff,              //!< effective plastic strian
        kelvin,             //!< absolute temperature
        DMG                 //!< cumulative damage for failure
    };
}

// Some inline function definition

inline void MPM3D_ErrorMessage(string filename, MPM_STATS linenumber, string msg)
{
    cout << "(MPM3D)Runtime Error in " << filename << " line " << linenumber << ":\n    " << msg << endl;
}

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