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
    Info: Implementation of class 'Failure_PriStress'
    Code-writter: Ruichen Ni
    Date: 2022.9.29
==============================================================*/

#include "Failure_PriStress.h"

Failure_PriStress::Failure_PriStress()
{
    Type = "Principle Stress";
    _min_principle_stress = 0.0;
    _max_principle_stress = 0.0;
    _max_shear_stress = 0.0;

    ParameterMap_Failure["PriStressMin"] = &_min_principle_stress;
    ParameterMap_Failure["PriStressMax"] = &_max_principle_stress;
    ParameterMap_Failure["ShearStressMax"] = &_max_shear_stress;
}

Failure_PriStress::~Failure_PriStress()
{
}

bool Failure_PriStress::CheckFailure(PhysicalProperty* pp)
{
    Array3D principle_stress;
    MPM_FLOAT mean_stress = pp->GetMeanStress();
    MPM_FLOAT bulk_q = pp->GetBulkViscosity();
    MPM_FLOAT stress_x = (*pp)[MPM::SDxx] + mean_stress - bulk_q;
    MPM_FLOAT stress_y = (*pp)[MPM::SDyy] + mean_stress - bulk_q;
    MPM_FLOAT stress_z = (*pp)[MPM::SDzz] + mean_stress - bulk_q;

    MPM_FLOAT I_1 = stress_x + stress_y + stress_z;
    MPM_FLOAT I_2 = -0.5*(I_1*I_1 - (stress_x*stress_x + stress_y*stress_y + stress_z*stress_z
        + 2.0*((*pp)[MPM::SDxy]*(*pp)[MPM::SDxy] + (*pp)[MPM::SDxz]*(*pp)[MPM::SDxz] +
               (*pp)[MPM::SDyz]*(*pp)[MPM::SDyz])));
    MPM_FLOAT I_3 = stress_x*stress_y*stress_z + 2.0*(*pp)[MPM::SDxy]*(*pp)[MPM::SDyz]*(*pp)[MPM::SDxz]
        - stress_x*(*pp)[MPM::SDyz]*(*pp)[MPM::SDyz] - stress_y*(*pp)[MPM::SDxz]*(*pp)[MPM::SDxz]
        - stress_z*(*pp)[MPM::SDxy]*(*pp)[MPM::SDxy];
    
    CubicFunctionRoots(1.0, -I_1, -I_2, -I_3, principle_stress);
    bool failure = false;
    MPM_FLOAT min_principle_stress = *min_element(principle_stress.begin(), principle_stress.end());
    MPM_FLOAT max_principle_stress = *max_element(principle_stress.begin(), principle_stress.end());
    MPM_FLOAT max_shear_stress = (max_principle_stress - min_principle_stress)*0.5;
    if (_min_principle_stress < -MPM_EPSILON &&
        min_principle_stress < -MPM_EPSILON &&
        min_principle_stress < _min_principle_stress)
        failure = true;
    
    if (_max_principle_stress > MPM_EPSILON &&
        max_principle_stress > MPM_EPSILON &&
        max_principle_stress > _max_principle_stress)
        failure = true;

    if (_max_shear_stress > MPM_EPSILON &&
        max_shear_stress > MPM_EPSILON &&
        max_shear_stress > _max_shear_stress)
        failure = true;
    
    if (failure)
    {
        pp->Failed();
        if(Erosion)
            pp->Eroded();
    }
    return pp->is_Failed();
}

void Failure_PriStress::Write(ofstream &os)
{
    os << "Failure model: " << Type << endl;
    os << "PriStressMin: " << _min_principle_stress << endl;
    os << "PriStressMax: " << _max_principle_stress << endl;
    os << "ShearStressMax: " << _max_shear_stress << endl;
    os << "Erosion: " << Erosion << endl;
}

bool Failure_PriStress::Initialize(map<string, MPM_FLOAT> &failure_para)
{
    if (!Failure_Base::Initialize(failure_para))
        return false;
    
    if (_min_principle_stress >= -MPM_EPSILON)
    {
        MPM3D_ErrorMessage(__FILE__, __LINE__, 
            "*** INPUT ERROR *** min_principle_strain should be less than zero.");
        return false;
    }

    if (_max_principle_stress <= MPM_EPSILON)
    {
        MPM3D_ErrorMessage(__FILE__, __LINE__, 
            "*** INPUT ERROR *** max_principle_strain should be greater than zero.");
        return false;
    }

    if (_max_shear_stress <= MPM_EPSILON)
    {
        MPM3D_ErrorMessage(__FILE__, __LINE__, 
            "*** INPUT ERROR *** max_shear_strain should be greater than zero.");
        return false;
    }
    return true;
}

bool Failure_PriStress::AddExtraParticleProperty_Failure(vector<MPM::ExtraParticleProperty> &ExtraProp)
{
    ExtraProp.push_back(MPM::SDxx);
    ExtraProp.push_back(MPM::SDxy);
    ExtraProp.push_back(MPM::SDxz);
    ExtraProp.push_back(MPM::SDyy);
    ExtraProp.push_back(MPM::SDyz);
    ExtraProp.push_back(MPM::SDzz);
    return true;
}