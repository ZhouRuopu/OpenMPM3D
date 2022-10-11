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
    Info: Class definition for Simple Johnson-Cook Plastic 
        strength behavior of material(LS-DYNA Material Model 98)
    Code-writter: Ruichen Ni
    Date: 2022.10.11
==============================================================*/

#ifndef _STRENGTH_SIMPLEJOHNSONCOOK_H_
#define _STRENGTH_SIMPLEJOHNSONCOOK_H_

#include "Strength_ElaPlastic.h"

//!> Simplified Johnson-Cook Model: No damage accumulated and no temperature in yield stress
class Strength_SimpleJohnsonCook: public Strength_ElaPlastic
{
public:
    Strength_SimpleJohnsonCook();
    ~Strength_SimpleJohnsonCook();

    //!> Write strength model information into file
    virtual void Write(ofstream &os);

    //!> Update the deviatoric stress of the particle
    virtual void UpdateDeviatoricStress(PhysicalProperty* pp, SymTensor& delta_strain, 
        SymTensor& delta_vortex, map<string, MPM_FLOAT>& transfer);
private:
    //!> Parameters for Johnson-Cook material
    MPM_FLOAT _B_jc, _n_jc, _C_jc;
    MPM_FLOAT _epso;    //!< strain rate normalization factor used in J-C model
};

#endif