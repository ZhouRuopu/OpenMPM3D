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
    Info: Implementation of class 'Strength_JohnsonCook'
    Code-writter: Ruichen Ni
    Date: 2022.10.11
==============================================================*/

#include "Strength_JohnsonCook.h"
#include "../../solver/Solver_Base.h"

Strength_JohnsonCook::Strength_JohnsonCook(/* args */)
{
    Type = "ISO-Plasticity: Johnson-Cook plasticity";
    _B_jc = 0.0;
    _n_jc = 0.0;
    _C_jc = 0.0;
    _m_jc = 0.0;
    _epso = 1.0;    //!< 1.0 for SI of Units
    _melt_temperature = 0.0;

    _compute_temperature = true;

    ParameterMap_Strength["B"] = &_B_jc;
    ParameterMap_Strength["n"] = &_n_jc;
    ParameterMap_Strength["C"] = &_C_jc;
    ParameterMap_Strength["m"] = &_m_jc;
    ParameterMap_Strength["melt"] = &_melt_temperature;
    ParameterMap_Strength["epso"] = &_epso;
}

Strength_JohnsonCook::~Strength_JohnsonCook()
{
}

void Strength_JohnsonCook::Write(ofstream& os)
{
    os << "Strength Type: " << Type << endl;
    os << "Density     E          Poisson" << endl;
    os << _density_0 << " " << _Young_Modulus << " " << _Poisson_Rate << endl;
    os << "Sigma_y       B          C          n          m" << endl;
    os << _yield_0 << " " << _B_jc << " " << _C_jc << " " << _n_jc << " " << _m_jc << endl;
    os << "Melt Temp.    Spec. Heat.    epso" << endl;
    os << _melt_temperature << " " << _specific_heat << " " << _epso << endl;
    if (_compute_temperature)
    {
        os << "Room Temp.    Temp. Coef.    Plas. Work Coef." << endl;
        os << _room_temperature << " " << _temperature_coefficient << " " << _plastic_work_coefficient << endl;
    }
    os << endl;
}

void Strength_JohnsonCook::UpdateDeviatoricStress(PhysicalProperty* pp, SymTensor& delta_strain, 
        SymTensor& delta_vortex, map<string, MPM_FLOAT>& transfer)
{
    _DamageDeviatoricStress(pp, delta_strain);
    pp->EquivalentStress();
    MPM_FLOAT seqv = pp->GetEquivalentStress();

    MPM_FLOAT depeff = 0.0;
    MPM_FLOAT lsrate = 0.0;
    MPM_FLOAT tstar = 0.0;
    transfer["yield"] = -1.0;
    MPM_FLOAT dt = Solver_Base::GetDTn_I();
    if (seqv > (*pp)[MPM::sigma_y])
    {
        tstar = ((*pp)[MPM::kelvin] - _room_temperature)/(_melt_temperature - _room_temperature);
        if (tstar > 1.0)    //!< melting
        {
            tstar = 1.0;
            pp->Failed();
        }
        else if (tstar < -MPM_EPSILON)
            tstar = 0.0;

        (*pp)[MPM::epeff] += 0.0001;    //!< avoid zero
        MPM_FLOAT plastic_modulus = _B_jc*_n_jc*pow((*pp)[MPM::epeff], _n_jc-1);    //!< Simplified
        (*pp)[MPM::epeff] -= 0.0001;
        depeff = (seqv - (*pp)[MPM::sigma_y])/(3.0*_shear_modulus + plastic_modulus);
        (*pp)[MPM::epeff] += depeff;

        MPM_FLOAT srate = depeff/_epso/dt;
        if (srate < 1.0)
            srate = 1.0;
        lsrate = log(srate);
        (*pp)[MPM::sigma_y] = (_yield_0 + _B_jc*pow((*pp)[MPM::epeff], _n_jc))*
            (1 + _C_jc*lsrate)*(1 - pow(tstar, _m_jc));
        if ((*pp)[MPM::sigma_y] < -MPM_EPSILON)
        {
            cout << "*** Warning: yield stress less than ZERO in Johnson-Cook Strength: "
                << (*pp)[MPM::sigma_y] << " with temperature of " << (*pp)[MPM::kelvin] << " K."<< endl;
        }

        if ((*pp)[MPM::sigma_y] > seqv)
        {
            (*pp)[MPM::epeff] -= depeff;
            depeff = 0.0;
        }
        else
        {
            transfer["yield"] = 1.0;
            MPM_FLOAT ratio = (*pp)[MPM::sigma_y]/seqv;
            pp->DeviatoricStressMultiplyScalar(ratio);
            pp->SetEquivalentStress((*pp)[MPM::sigma_y]);
        }
    }
    transfer["depeff"] = depeff;
    transfer["lsrate"] = lsrate;
    transfer["tstar"] = tstar;
}

bool Strength_JohnsonCook::AddExtraParticleProperty_Strength(vector<MPM::ExtraParticleProperty> &ExtraProp,
        map<string, MPM_FLOAT>& transfer)
{
    if (!Strength_ElaPlastic::AddExtraParticleProperty_Strength(ExtraProp, transfer))
        return false;
    ExtraProp.push_back(MPM::DMG);
    return true;
}