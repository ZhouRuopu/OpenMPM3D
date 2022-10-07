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
    Info: Base class definition for all solvers
    Code-writter: Ruichen Ni
    Date: 2022.9.28
==============================================================*/

#ifndef _SOLVER_BASE_H
#define _SOLVER_BASE_H

#include "../main/MPM3D_MACRO.h"

class Solver_Base
{
public:
    Solver_Base();
    ~Solver_Base();
protected:
    //!> Global variables which can be obtained by static Get() function
    static MPM_FLOAT _dtn,             //!< Time step: t^(n-1/2) = t^n - t^(n-1)
                     _dtn1,            //!< Time step: t^(n+1/2) = t^(n+1) - t^n
                     _dtn1_half,       //!< Time step: _dtn1*0.5
                     _dtx,             //!< (_dtn + _dtn1)*0.5
                     _current_time;

//!> Various Get/Set function
public:
    inline static MPM_FLOAT GetDTn() {return _dtn;}
    inline static MPM_FLOAT GetDTn_I() {return _dtn1;}
    inline static MPM_FLOAT GetDTn_I_Half() {return _dtn1_half;}
    inline static MPM_FLOAT GetDTx() {return _dtx;}
    inline static MPM_FLOAT GetCurrentTime() {return _current_time;}
};

#endif