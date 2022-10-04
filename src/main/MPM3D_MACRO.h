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
#include <map>
#include <iostream>
#include <fstream>
using namespace std;

#define __min(x, y) (x < y?x:y)
#define __max(x, y) (x > y?x:y)
#define __min3(x, y, z) (__min(__min(x, y), z))
#define __max3(x, y, z) (__max(__max(x, y), z))

#define MPM3D_ErrorMessage(msg) \
    cout << "(MPM3D)Runtime Error in " << __FILE__ << " line " << __LINE__ << ":\n    " << msg << "\n"

#ifdef _MPM_DOUBLE
    #define MPM_FLOAT       double
    #define MPM_EPSILON     DBL_EPSILON
    #define MPM_MIN         DBL_MIN
    #define MPM_MAX         DBL_MAX
    #define PI              3.14159265358979323846
    #define SQRT3           1.7320508075688773
#else
    #define MPM_FLOAT       float
    #define MPM_EPSILON     FLT_EPSILON
    #define MPM_MIN         FLT_MIN
    #define MPM_MAX         FLT_MAX
    #define PI              3.1415926536F
    #define SQRT3           1.7320508076F
#endif

#ifdef _MPM_MASSIVE_PARTICLE
    #define MPM_STATS       unsigned long
#else
    #define MPM_STATS       unsigned int
#endif

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
    const int ExtraParticlePropertySum = 15;
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
#endif