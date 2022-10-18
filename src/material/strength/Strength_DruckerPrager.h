#pragma once
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
    Info: Class definition for Dracker-Prager strength
        behavior of material(LS-DYNA Material Model 193)
    Code-writter: Ruopu Zhou
    Date: 2022.10.17
==============================================================*/

#ifndef _STRENGTH_DRUCKERPRAGER_H_
#define _STRENGTH_DRUCKERPRAGER_H_

#include "Strength_ElaPlastic.h"

class Strength_DruckerPrager : public Strength_ElaPlastic
{
public:
    Strength_DruckerPrager();
    ~Strength_DruckerPrager();

    //!> Write strength model information into file
    virtual void Write(ofstream& os);

    //!> Update the deviatoric stress of the particle
    virtual void UpdateDeviatoricStress(PhysicalProperty* pp, SymTensor& delta_strain,
        SymTensor& delta_vortex, map<string, MPM_FLOAT>& transfer);

    //!> Update the pressure of the particle with elastic assumption
    virtual void ElasticPressure(PhysicalProperty* pp, MPM_FLOAT delta_vol,
        map<string, MPM_FLOAT>& transfer);

private:
    //!> Parameters for Johnson-Cook material
    MPM_FLOAT _q_fai, _k_fai, _q_psi, _ten_f;

protected:
    MPM_FLOAT iplas;
};

#endif