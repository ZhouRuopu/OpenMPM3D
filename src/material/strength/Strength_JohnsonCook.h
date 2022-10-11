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
    Info: Class definition for Johnson-Cook Plastic strength 
        behavior of material(LS-DYNA Material Model 15)
    Code-writter: Ruichen Ni
    Date: 2022.10.11
==============================================================*/

#ifndef _STRENGTH_JOHNSONCOOK_H_
#define _STRENGTH_JOHNSONCOOK_H_

#include "Strength_ElaPlastic.h"

class Strength_JohnsonCook: public Strength_ElaPlastic
{
public:
    Strength_JohnsonCook();
    ~Strength_JohnsonCook();

    //!> Write strength model information into file
    virtual void Write(ofstream &os);

    //!> Update the deviatoric stress of the particle
    virtual void UpdateDeviatoricStress(PhysicalProperty* pp, SymTensor& delta_strain, 
        SymTensor& delta_vortex, map<string, MPM_FLOAT>& transfer);

    //!> Add extra particle properties based on different strength model
    virtual bool AddExtraParticleProperty_Strength(vector<MPM::ExtraParticleProperty> &ExtraProp,
        map<string, MPM_FLOAT>& transfer);
private:
    //!> Parameters for Johnson-Cook material
    MPM_FLOAT _B_jc, _n_jc, _C_jc, _m_jc;
    MPM_FLOAT _epso;    //!< strain rate normalization factor used in J-C model
    MPM_FLOAT _melt_temperature;
};

#endif