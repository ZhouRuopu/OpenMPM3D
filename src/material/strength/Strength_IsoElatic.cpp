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
    MPM_FLOAT delta_strain_mean = (delta_strain[0] + delta_strain[1] + delta_strain[2])/3.0;

    SymTensor deviatoric_stress = pp->GetDeviatoricStress();
    deviatoric_stress[0] += 2.0*_shear_modulus*(delta_strain[0] - delta_strain_mean);
    deviatoric_stress[1] += 2.0*_shear_modulus*(delta_strain[1] - delta_strain_mean);
    deviatoric_stress[2] += 2.0*_shear_modulus*(delta_strain[2] - delta_strain_mean);
    deviatoric_stress[3] += _shear_modulus*delta_strain[3];
    deviatoric_stress[4] += _shear_modulus*delta_strain[4];
    deviatoric_stress[5] += _shear_modulus*delta_strain[5];
    pp->SetDeviatoricStress(deviatoric_stress);

    pp->EquivalentStress();
}

void Strength_IsoElastic::ElasticPressure(PhysicalProperty* pp, MPM_FLOAT delta_vol,
        map<string, MPM_FLOAT>& transfer)
{
    if (_compute_temperature)
    {
        MPM_FLOAT origin_volume = pp->GetMass()/_density_0;
        MPM_FLOAT mean_stress = _volumetric_modulus*(pp->GetVolume() - origin_volume)/origin_volume - 
            _temperature_coefficient*((*pp)[MPM::kelvin] - _room_temperature);
        pp->SetMeanStress(mean_stress);
    }
    else
    {
        MPM_FLOAT mean_stress = pp->GetMeanStress();
        mean_stress += _volumetric_modulus*delta_vol;
        pp->SetMeanStress(mean_stress);
    }
}

void Strength_IsoElastic::UpdateTemperature(PhysicalProperty* pp, bool yield, MPM_FLOAT delta_vol,
        map<string, MPM_FLOAT>& transfer)
{
    if (_compute_temperature)
        (*pp)[MPM::kelvin] -= _temperature_coefficient*(*pp)[MPM::kelvin]*delta_vol/pp->GetDensity()/_specific_heat;
}

bool Strength_IsoElastic::AddExtraParticleProperty_Strength(vector<MPM::ExtraParticleProperty> &ExtraProp)
{
    if (_compute_temperature)
        ExtraProp.push_back(MPM::kelvin);
    return true;
}