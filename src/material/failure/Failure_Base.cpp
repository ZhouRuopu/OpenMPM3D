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
    Info: Implementation of class 'Failure_Base'
    Code-writter: Ruichen Ni
    Date: 2022.9.28
==============================================================*/

#include "./Failure_Base.h"

void CubicFunctionSolution(MPM_FLOAT a, MPM_FLOAT b, MPM_FLOAT c, MPM_FLOAT d, MPM_FLOAT (&roots)[3])
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

Failure_Base::Failure_Base()
{
    Type = "";
    Erosion = false;
}

Failure_Base::~Failure_Base()
{
}

bool Failure_Base::Initialize(map<string, MPM_FLOAT> &failure_para)
{
    for(map<string, MPM_FLOAT>::iterator iter = failure_para.begin(); 
        iter != failure_para.end(); iter++)
    {
        if(ParameterMap_Failure.find(iter->first) != ParameterMap_Failure.end())
            *ParameterMap_Failure[iter->first] = iter->second;
        else
        {
            string error_msg = "Can't find the failure model parameter " + iter->first + " at " + Type;
            MPM3D_ErrorMessage(error_msg);
            return false;
        }
    }
    return true;
}