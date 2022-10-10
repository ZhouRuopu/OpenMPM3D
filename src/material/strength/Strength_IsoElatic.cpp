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
    Info: Implementation of class 'Strength_IsoElastic'
    Code-writter: Ruichen Ni
    Date: 2022.10.10
==============================================================*/

#include "Strength_IsoElastic.h"

Strength_IsoElastic::Strength_IsoElastic()
{
    Type = "ISO-Elasticity: Isotropic elasticity";
}

Strength_IsoElastic::~Strength_IsoElastic()
{
}

void Strength_IsoElastic::Write(ofstream& os)
{
    os << "Strength Type: " << Type << endl;
    os << "Density     E          Poisson" << endl;
    os << _density_0 << " " << _Young_Modulus << " " << _Poisson_Rate << endl;
    if (_compute_temperature)
    {
        os << "Room Temp.    Temp. Coef." << endl;
        os << _room_temperature << " " << _temperature_coefficient << endl;
    }
    os << endl;
}

void Strength_IsoElastic::UpdateDeviatoricStress(PhysicalProperty* pp, SymTensor& delta_strain, 
        SymTensor& delta_vortex, map<string, MPM_FLOAT>& transfer)
{
    _ElasticDeviatoricStress(pp, delta_strain);

    pp->EquivalentStress();
}

void Strength_IsoElastic::ElasticPressure(PhysicalProperty* pp, MPM_FLOAT delta_vol,
        map<string, MPM_FLOAT>& transfer)
{
    _ElasticPressure(pp, delta_vol);
}

void Strength_IsoElastic::UpdateTemperature(PhysicalProperty* pp, MPM_FLOAT delta_vol,
        map<string, MPM_FLOAT>& transfer)
{
    if (_compute_temperature)
        (*pp)[MPM::kelvin] -= _temperature_coefficient*(*pp)[MPM::kelvin]*delta_vol/pp->GetDensity()/_specific_heat;
}