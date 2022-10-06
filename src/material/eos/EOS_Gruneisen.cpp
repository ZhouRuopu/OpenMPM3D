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
    Info: Implementation of class "EOS_Gruneisen"
    Code-writter: Ruichen Ni
    Date: 2022.10.5
==============================================================*/

#include "EOS_Gruneisen.h"

EOS_Gruneisen::EOS_Gruneisen()
{
    Type = "Mie-Gruneisen EOS";
    _cutoff = 0.9;

    ParameterMap_EOS["S1"] = &_s;
    ParameterMap_EOS["gamma0"] = &_gamma0;
    ParameterMap_EOS["cutoff"] = &_cutoff;
}

EOS_Gruneisen::~EOS_Gruneisen()
{
}

bool EOS_Gruneisen::Initialize(map<string, MPM_FLOAT> &eos_para, MPM_FLOAT rho0)
{
    if (!EOS_Base::Initialize(eos_para, rho0))
        return false;
    
    if (_sound_speed_0 < MPM_EPSILON)
    {
        string error_msg = "*** Error *** Mie-Gruneisen EOS parameter C0 is needed !";
        MPM3D_ErrorMessage(__FILE__, __LINE__, error_msg);
        return false;
    }

    _impendence_0 = _density_0*_sound_speed_0*_sound_speed_0;
    _mu_cutoff = _cutoff/(_s - 1.0);

    MPM_FLOAT denominator = 1.0 - (_s - 1.0)*_mu_cutoff;
    _pH_cutoff = _impendence_0*_mu_cutoff*(1.0 + _mu_cutoff)/(denominator*denominator);
    _DpH_cutoff = _impendence_0*(1.0 + (_s + 1.0)*_mu_cutoff)/(denominator*denominator*denominator);
    
    return true;
}

void EOS_Gruneisen::Write(ofstream& os)
{
    os << "EOS Type: " << Type << endl;
    os << "C0          S1          Gamma0          Cutoff" << endl;
    os << _sound_speed_0 << " " << _s << " " << _gamma0 << " " << _cutoff << endl << endl;
}

void EOS_Gruneisen::UpdatePressure(PhysicalProperty* pp, MPM_FLOAT delta_vol_half, 
    MPM_FLOAT delta_ie, map<string, MPM_FLOAT>& transfer)
{
    MPM_FLOAT V0 = pp->GetMass()/_density_0;
    MPM_FLOAT E = (pp->GetInternalEnergy() + delta_ie)/V0;
    MPM_FLOAT mu = pp->GetDensity()/_density_0 - 1.0;
    MPM_FLOAT gamma = _gamma0*_density_0/pp->GetDensity();

    MPM_FLOAT A;
    if (mu > MPM_EPSILON)
    {
        MPM_FLOAT pH;
        if (mu > _mu_cutoff)
        {
            pH = _pH_cutoff + _DpH_cutoff*(mu - _mu_cutoff);
            pp->Failed();
        }
        else
        {
            MPM_FLOAT denominator = 1.0 - (_s - 1.0)*mu;
            pH = _impendence_0*mu*(1.0 + mu)/(denominator*denominator);
        }
        A = pH*(1 - 0.5*gamma*mu);
    }
    else
    {
        A = _impendence_0*mu;
    }

    MPM_FLOAT B = _gamma0;

    MPM_FLOAT pressure_new = (A + B*E)/(1 + B*delta_vol_half/V0);
    pp->SetMeanStress(-pressure_new);
    return;
}

MPM_FLOAT EOS_Gruneisen::SoundSpeedSquare_EOS(PhysicalProperty* pp)
{
    if (pp->is_Failed())
        return 0.0;
    
    MPM_FLOAT pressure = -pp->GetMeanStress();
    MPM_FLOAT rv = _density_0/pp->GetDensity();     //!< Relative volume
    MPM_FLOAT mu = 1.0/rv - 1.0;
    MPM_FLOAT gamma = _gamma0*rv;

    MPM_FLOAT result;
    if (mu > MPM_EPSILON)
    {
        MPM_FLOAT pH, DpH;
        MPM_FLOAT denominator = 1.0 - (_s - 1.0)*mu;
        pH = _impendence_0*mu*(1.0 + mu)/(denominator*denominator);
        DpH = _impendence_0*(1.0 + (_s + 1.0)*mu)/(denominator*denominator*denominator);

        result = (DpH*(1 - 0.5*gamma*mu) - 0.5*pH*gamma)/_density_0 
            + pressure*rv*gamma/_density_0;
    }
    else
    {
        result = _sound_speed_0*_sound_speed_0 + pressure*rv*gamma/_density_0;
    }
    return result;
}