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
    Info: Implementation of class 'Strength_SimpleJohnsonCook'
    Code-writter: Ruichen Ni
    Date: 2022.10.11
==============================================================*/

#include "Strength_SimpleJohnsonCook.h"
#include "../../solver/Solver_Base.h"

Strength_SimpleJohnsonCook::Strength_SimpleJohnsonCook()
{
    Type = "ISO-Plasticity: Simplified Johnson-Cook plasticity";
    _B_jc = 0.0;
    _n_jc = 0.0;
    _C_jc = 0.0;
    _epso = 1.0;    //!< 1.0 for SI of Units

    ParameterMap_Strength["B"] = &_B_jc;
    ParameterMap_Strength["n"] = &_n_jc;
    ParameterMap_Strength["C"] = &_C_jc;
    ParameterMap_Strength["epso"] = &_epso;
}

Strength_SimpleJohnsonCook::~Strength_SimpleJohnsonCook()
{
}

void Strength_SimpleJohnsonCook::Write(ofstream& os)
{
    os << "Strength Type: " << Type << endl;
    os << "Density     E          Poisson" << endl;
    os << _density_0 << " " << _Young_Modulus << " " << _Poisson_Rate << endl;
    os << "Sigma_y       B          C          n          epso" << endl;
    os << _yield_0 << " " << _B_jc << " " << _C_jc << " " << _n_jc << " " << _epso << endl;
    if (_compute_temperature)
    {
        os << "Room Temp.    Temp. Coef.    Plas. Work Coef." << endl;
        os << _room_temperature << " " << _temperature_coefficient << " " << _plastic_work_coefficient << endl;
    }
    os << endl;
}

void Strength_SimpleJohnsonCook::UpdateDeviatoricStress(PhysicalProperty* pp, SymTensor& delta_strain, 
        SymTensor& delta_vortex, map<string, MPM_FLOAT>& transfer)
{
    _ElasticDeviatoricStress(pp, delta_strain);
    pp->EquivalentStress();
    MPM_FLOAT seqv = pp->GetEquivalentStress();

    MPM_FLOAT depeff = 0.0;
    transfer["yield"] = -1.0;
    MPM_FLOAT dt = Solver_Base::GetDTn_I();
    if (seqv > (*pp)[MPM::sigma_y])
    {
        (*pp)[MPM::epeff] += 0.0001;    //!< avoid zero
        MPM_FLOAT plastic_modulus = _B_jc*_n_jc*pow((*pp)[MPM::epeff], _n_jc-1);    //!< Simplified
        (*pp)[MPM::epeff] -= 0.0001;
        depeff = (seqv - (*pp)[MPM::sigma_y])/(3.0*_shear_modulus + plastic_modulus);
        (*pp)[MPM::epeff] += depeff;

        MPM_FLOAT srate = depeff/_epso/dt;
        if (srate < 1.0)
            srate = 1.0;
        (*pp)[MPM::sigma_y] = (_yield_0 + _B_jc*pow((*pp)[MPM::epeff], _n_jc))*
            (1 + _C_jc*log(srate));
        if ((*pp)[MPM::sigma_y] < -MPM_EPSILON)
        {
            cout << "*** Warning: yield stress less than ZERO in Simplified Johnson-Cook Strength: "
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
}