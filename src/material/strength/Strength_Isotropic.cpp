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
    Info: Implementation of class 'Strength_Isotropic'
    Code-writter: Ruichen Ni
    Date: 2022.10.9
==============================================================*/

#include "Strength_Isotropic.h"

Strength_Isotropic::Strength_Isotropic()
{
    Type = "";  //!< Strength_Isotropic cannot be Instantiated

    _Young_Modulus = 0.0;
    _Poisson_Rate = 0.0;
    _shear_modulus = 0.0;
    _volumetric_modulus = 0.0;

    _specific_heat = 0.0;
    _temperature_coefficient = 0.0;
    _room_temperature = 293.0;

    ParameterMap_Strength["Young"] = &_Young_Modulus;
    ParameterMap_Strength["Poisson"] = &_Poisson_Rate;
    ParameterMap_Strength["SpecHeat"] = &_specific_heat;
    ParameterMap_Strength["TemperatureCoefficient"] = &_temperature_coefficient;
    ParameterMap_Strength["roomt"] = &_room_temperature;
}

Strength_Isotropic::~Strength_Isotropic()
{
}

bool Strength_Isotropic::Initialize(map<string, MPM_FLOAT> &strength_para, MPM_FLOAT rho0)
{
    if (!Strength_Base::Initialize(strength_para, rho0))
        return false;
    
    _shear_modulus = 0.5*_Young_Modulus/(1.0 + _Poisson_Rate);
    _volumetric_modulus = _Young_Modulus/(3.0*(1.0 - 2.0*_Poisson_Rate));
    return true;
}

MPM_FLOAT Strength_Isotropic::SoundSpeedSquare_Strength(PhysicalProperty* pp)
{
    return FourThird*_shear_modulus/pp->GetDensity();
}

MPM_FLOAT Strength_Isotropic::SoundSpeedSquare_Elastic(PhysicalProperty* pp)
{
    return _volumetric_modulus/pp->GetDensity();
}