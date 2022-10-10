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
    Info: Implementation of class 'Strength_ElaPlastic'
    Code-writter: Ruichen Ni
    Date: 2022.10.10
==============================================================*/

#include "Strength_ElaPlastic.h"

Strength_ElaPlastic::Strength_ElaPlastic()
{
    Type = "ISO-Plasticity: Elastic-perfectly plasticity";

    _plastic_work_coefficient = 0.9;

    ParameterMap_Strength["Yield0"] = &_yield_0;
    ParameterMap_Strength["PlasticWorkCoefficient"] = &_plastic_work_coefficient;
}

Strength_ElaPlastic::~Strength_ElaPlastic()
{
}

void Strength_ElaPlastic::Write(ofstream& os)
{
    os << "Strength Type: " << Type << endl;
    os << "Density     E          Poisson        Sigma_y" << endl;
    os << _density_0 << " " << _Young_Modulus << " " << _Poisson_Rate << " " << _yield_0 << endl;
    if (_compute_temperature)
    {
        os << "Room Temp.    Temp. Coef." << endl;
        os << _room_temperature << " " << _temperature_coefficient << endl;
    }
    os << endl;
}

void Strength_ElaPlastic::UpdateDeviatoricStress(PhysicalProperty* pp, SymTensor& delta_strain, 
        SymTensor& delta_vortex, map<string, MPM_FLOAT>& transfer)
{
    _ElasticDeviatoricStress(pp, delta_strain);

    pp->EquivalentStress();
    MPM_FLOAT seqv = pp->GetEquivalentStress();

    MPM_FLOAT depeff = 0.0;
    bool yield = (seqv > _yield_0);
    transfer["yield"] = -1.0;   //!< means not yield, used for temperature update
    if (yield)
    {
        transfer["yield"] = 1.0;

        depeff = (seqv - (*pp)[MPM::sigma_y])/(3.0*_shear_modulus);
        transfer["depeff"] = depeff;

        (*pp)[MPM::epeff] += depeff;
        MPM_FLOAT ratio = _yield_0/seqv;
        pp->DeviatoricStressMultiplyScalar(ratio);
        pp->SetEquivalentStress(seqv*ratio);
    }
}

void Strength_ElaPlastic::ElasticPressure(PhysicalProperty* pp, MPM_FLOAT delta_vol,
        map<string, MPM_FLOAT>& transfer)
{
    _ElasticPressure(pp, delta_vol);
}

void Strength_ElaPlastic::UpdateTemperature(PhysicalProperty* pp, MPM_FLOAT delta_vol,
        map<string, MPM_FLOAT>& transfer)
{
    if (_compute_temperature)
    {
        if (pp->is_Failed())
            return;
        
        if (transfer["yield"] < 0.0)    //!< must be available because of  "UpdateDeviatoricStress" function
            return;
        
        (*pp)[MPM::kelvin] += _plastic_work_coefficient*pp->GetEquivalentStress()*
            transfer["depeff"]/pp->GetDensity()/_specific_heat;
    }
}

bool Strength_ElaPlastic::AddExtraParticleProperty_Strength(vector<MPM::ExtraParticleProperty> &ExtraProp,
        map<string, MPM_FLOAT>& transfer)
{
    if (!Strength_Isotropic::AddExtraParticleProperty_Strength(ExtraProp, transfer))
        return false;
    ExtraProp.push_back(MPM::epeff);
    ExtraProp.push_back(MPM::sigma_y);
    transfer["sigma_y"] = _yield_0;     //!< for particle property "sigma_y" initialization
    return true;
}