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
    Info: Implementation of class "Solver_Base"
    Code-writter: Ruichen Ni
    Date: 2022.9.28
==============================================================*/

#include "Solver_Base.h"

MPM_FLOAT Solver_Base::_dtn = 0.0;
MPM_FLOAT Solver_Base::_dtn1 = 0.0;
MPM_FLOAT Solver_Base::_dtn1_half = 0.0;
MPM_FLOAT Solver_Base::_current_time = 0.0;

Solver_Base::Solver_Base()
{
}

Solver_Base::~Solver_Base()
{
}