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
    Info: Base class definition for Isotropic strength behavior 
        of material
    Code-writter: Ruichen Ni
    Date: 2022.10.9
==============================================================*/

#ifndef _STRENGTH_ISOTROPIC_H_
#define _STRENGTH_ISOTROPIC_H_

#include "Strength_Base.h"

class Strength_Isotropic:public Strength_Base
{
public:
    Strength_Isotropic();
    ~Strength_Isotropic();

    //!> Initial the strength model with parameters' map
    virtual bool Initialize(map<string, MPM_FLOAT> &strength_para, MPM_FLOAT rho0);

    //!> Write strength model information into file
    virtual void Write(ofstream &os) = 0;

    //!> Update the deviatoric stress of the particle
    //!> the map named transfer is used to tansfer data between Strength/EOS/Failure model
    virtual void UpdateDeviatoricStress(PhysicalProperty* pp, SymTensor& delta_strain, 
        SymTensor& delta_vortex, map<string, MPM_FLOAT>& transfer) = 0;

    //!> Update the pressure of the particle with elastic assumption
    virtual void ElasticPressure(PhysicalProperty* pp, MPM_FLOAT delta_vol,
        map<string, MPM_FLOAT>& transfer) = 0;

    //!> Calculate the squared adabatic sound speed of deviatoric part
    virtual MPM_FLOAT SoundSpeedSquare_Strength(PhysicalProperty* pp);

    //!> Calculate the squared adabatic sound speed of elastic part
    virtual MPM_FLOAT SoundSpeedSquare_Elastic(PhysicalProperty* pp);

    //!> Update the temperature of the particle
    virtual void UpdateTemperature(PhysicalProperty* pp, bool yield, MPM_FLOAT delta_vol,
        map<string, MPM_FLOAT>& transfer) = 0;

    //!> Add extra particle properties based on different strength model
    virtual bool AddExtraParticleProperty_Strength(vector<MPM::ExtraParticleProperty> &ExtraProp) = 0;
protected:
    //!> Isotropic constitution consists of two parameters
    //!>    Young's Modulus and Poisson Rate
    MPM_FLOAT _Young_Modulus;
    MPM_FLOAT _Poisson_Rate;

    //!> Other parameters obtained by calculations of Young's Modulus and Poisson Rate
    MPM_FLOAT _shear_modulus;
    MPM_FLOAT _volumetric_modulus;

    //!> Parameters corresponding temperature
    MPM_FLOAT _specific_heat;
    MPM_FLOAT _temperature_coefficient;
    MPM_FLOAT _room_temperature;
};

#endif