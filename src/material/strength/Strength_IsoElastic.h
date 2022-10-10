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
    Info: Class definition for Isotropic Elastic strength 
        behavior of material
    Code-writter: Ruichen Ni
    Date: 2022.10.10
==============================================================*/

#ifndef _STRENGTH_ISOELASTIC_H_
#define _STRENGTH_ISOELASTIC_H_

#include "Strength_Isotropic.h"

class Strength_IsoElastic: public Strength_Isotropic
{
public:
    Strength_IsoElastic();
    ~Strength_IsoElastic();

    //!> Write strength model information into file
    virtual void Write(ofstream &os);

    //!> Update the deviatoric stress of the particle
    //!> the map named transfer is used to tansfer data between Strength/EOS/Failure model
    virtual void UpdateDeviatoricStress(PhysicalProperty* pp, SymTensor& delta_strain, 
        SymTensor& delta_vortex, map<string, MPM_FLOAT>& transfer);

    //!> Update the pressure of the particle with elastic assumption
    virtual void ElasticPressure(PhysicalProperty* pp, MPM_FLOAT delta_vol,
        map<string, MPM_FLOAT>& transfer);

    //!> Update the temperature of the particle
    virtual void UpdateTemperature(PhysicalProperty* pp, bool yield, MPM_FLOAT delta_vol,
        map<string, MPM_FLOAT>& transfer);

    //!> Add extra particle properties based on different strength model
    virtual bool AddExtraParticleProperty_Strength(vector<MPM::ExtraParticleProperty> &ExtraProp);
protected:
    /* data */
};

#endif