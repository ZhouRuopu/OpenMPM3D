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
    Info: Implementation of class 'Strength_Drucker-Prager'
    Code-writter: Ruopu Zhou
    Date: 2022.10.17
==============================================================*/

#include "Strength_DruckerPrager.h"
#include "../../solver/Solver_Base.h"

Strength_DruckerPrager::Strength_DruckerPrager()
{
    Type = "ISO-Plasticity: DruckerPrager plasticity";
    _q_fai = 0.0;
    _k_fai = 0.0;
    _q_psi = 0.0;
    ParameterMap_Strength["qfai"] = &_q_fai;
    ParameterMap_Strength["kfai"] = &_k_fai;
    ParameterMap_Strength["qpsi"] = &_q_psi;
    ParameterMap_Strength["tenf"] = &_ten_f;
}

Strength_DruckerPrager::~Strength_DruckerPrager()
{
}

void Strength_DruckerPrager::Write(ofstream& os)
{
    os << "Strength Type: " << Type << endl;
    os << "Density     E          Poisson" << endl;
    os << _density_0 << " " << _Young_Modulus << " " << _Poisson_Rate << endl;
    os << "qfai          kfai          qpsi          tenf" << endl;
    os << _q_fai << " " << _k_fai << " " << _q_psi << " " << _ten_f << " " << endl;
    os << endl;
}

void Strength_DruckerPrager::UpdateDeviatoricStress(PhysicalProperty* pp, SymTensor& delta_strain,
    SymTensor& delta_vortex, map<string, MPM_FLOAT>& transfer)
{
    _ElasticDeviatoricStress(pp, delta_strain);
    pp->EquivalentStress();
    MPM_FLOAT seqv = pp->GetEquivalentStress(); //EquivalentStress
    MPM_FLOAT _new_mean_stress = pp->GetMeanStress(); //MeanStress

    MPM_FLOAT depeff = 0.0;
    MPM_FLOAT _G_modulus = _Young_Modulus / (2 * (1 + _Poisson_Rate));
    MPM_FLOAT _K_modulus = _Young_Modulus / (3 * (1 - _Poisson_Rate));
    MPM_FLOAT _Tau = seqv/sqrt(3.0);
    
    MPM_FLOAT tenf = 0.0; //Degenerate to Mises yield surface
    if (_q_fai != 0.0)
    {
        MPM_FLOAT tenf_max = _k_fai / _q_fai;
        tenf = min(_ten_f, tenf_max);
    }
    MPM_FLOAT _DP_Fi = _Tau + _q_fai * _new_mean_stress - _k_fai; //D-P yeild surface
    MPM_FLOAT _DP_sig = _new_mean_stress - tenf; // the spherical stress difference

    MPM_FLOAT _dlamd = 0.0;
    MPM_FLOAT _new_Tau = 0.0;
    MPM_FLOAT ratio = 0.0;
    MPM_FLOAT _DP_hfai = 0.0;
    MPM_FLOAT _Tau_p = 0.0;
    MPM_FLOAT _alpha_p = 0.0;
    transfer["yield"] = -1.0;
    if (_DP_sig <-MPM_EPSILON)
    {
        if (_DP_Fi > MPM_EPSILON)
        {
            iplas = 1;
            transfer["yield"] = 1.0;
            _dlamd = _DP_Fi / (_G_modulus + _K_modulus * _q_fai * _q_psi);
            _new_mean_stress = _new_mean_stress - _K_modulus * _q_psi * _dlamd;
            _new_Tau = _k_fai - _q_fai * _new_mean_stress;
            ratio = _new_Tau / _Tau;

            pp->DeviatoricStressMultiplyScalar(ratio);
            pp->SetEquivalentStress(seqv * ratio);
            pp->SetMeanStress(_new_mean_stress);

            depeff = _dlamd * sqrt(1.0 / 3.0 + (2.0 / 9.0) * pow(_q_psi, 2));
            transfer["depeff"] = depeff;
            (*pp)[MPM::epeff] += depeff;
        }
    }
    else //_DP_sig >= 0.0
    {
        _alpha_p = sqrt(1 + pow(_q_fai, 2)) - _q_fai;
        _Tau_p = _k_fai - _q_fai * tenf;
        _DP_hfai = _Tau - _Tau_p - _alpha_p * _DP_sig;

        if (_DP_hfai > MPM_EPSILON) // same as _DP_Fi > 0.0
        {
            iplas = 1;
            transfer["yield"] = 1.0;
            _dlamd = _DP_Fi / (_G_modulus + _K_modulus * _q_fai * _q_psi);
            _new_mean_stress = _new_mean_stress - _K_modulus * _q_psi * _dlamd;
            _new_Tau = _k_fai - _q_fai * _new_mean_stress;
            ratio = _new_Tau / _Tau;

            pp->DeviatoricStressMultiplyScalar(ratio);
            pp->SetEquivalentStress(seqv * ratio);
            pp->SetMeanStress(_new_mean_stress);

            depeff = _dlamd * sqrt(1.0 / 3.0 + (2.0 / 9.0) * pow(_q_psi, 2));
            transfer["depeff"] = depeff;
            (*pp)[MPM::epeff] += depeff;
        }
        else // _DP_hfai <= 0.0
        {
            iplas = 2;
            transfer["yield"] = 1.0;
            _dlamd = (_new_mean_stress - tenf) / _K_modulus;
            pp->SetMeanStress(tenf); // update Mean Stress only
            
            depeff = _dlamd * (1.0 / 3.0) * sqrt(2.0);
            transfer["depeff"] = depeff;
            (*pp)[MPM::epeff] += depeff;
        }
    }

}

void Strength_DruckerPrager::ElasticPressure(PhysicalProperty* pp, MPM_FLOAT delta_vol,
    map<string, MPM_FLOAT>& transfer)
{
    _ElasticPressure(pp, delta_vol);
}