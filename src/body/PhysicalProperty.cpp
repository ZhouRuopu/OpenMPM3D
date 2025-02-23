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
    Info: Implementation of class 'PhysicalProperty'
    Code-writter: Ruichen Ni
    Date: 2022.9.29
==============================================================*/

#include "./PhysicalProperty.h"
#include "../utility/MathFunctionList.h"

PhysicalProperty::PhysicalProperty()
{
    _mass = 0.0;
    _volume = 0.0;
    _density = 0.0;
    _mean_stress = 0.0;
    _deviatoric_stress.fill(0.0);
    _equivalent_stress = 0.0;
    _bulk_q = 0.0;
    _internal_energy = 0.0;
    _sound_speed = 0.0;
    _failure = false;
    _eroded = false;

    _extra_properties = nullptr;
    _extra_property_positions = nullptr;
}

PhysicalProperty::PhysicalProperty(const PhysicalProperty& pp)
{
    _mass = pp._mass;
    _volume = pp._volume;
    _density = pp._density;
    _mean_stress = pp._mean_stress;
    _deviatoric_stress = pp._deviatoric_stress;
    _equivalent_stress = pp._equivalent_stress;
    _bulk_q = pp._bulk_q;
    _internal_energy = pp._internal_energy;
    _sound_speed = pp._sound_speed;
    _failure = pp._failure;
    _eroded = pp._eroded;

    _extra_properties = nullptr;
    _extra_property_positions = pp._extra_property_positions;
    if (_extra_property_positions)
    {
        int count = 0;
        for (int i = 0; i < MPM::ExtraParticlePropertySum; i++)
            if (_extra_property_positions[i] > 0)
                count++;
        
        AllocateMemoryForExtraParticleProperty(count);
        for (int i = 0; i < count; i++)
            _extra_properties[i] = pp._extra_properties[i];
    }
}

PhysicalProperty::~PhysicalProperty()
{
    if (_extra_properties)
        delete[] _extra_properties;
}

void PhysicalProperty::AllocateMemoryForExtraParticleProperty(int number)
{
    _extra_properties = new MPM_FLOAT[number];
    for (int i = 0; i < number; i++)
        _extra_properties[i] = 0.0;
}

Array3D&& PhysicalProperty::CalculatePrincipleStress()
{
    MPM_FLOAT stress_x = _deviatoric_stress[0] + _mean_stress - _bulk_q;
    MPM_FLOAT stress_y = _deviatoric_stress[1] + _mean_stress - _bulk_q;
    MPM_FLOAT stress_z = _deviatoric_stress[2] + _mean_stress - _bulk_q;

    MPM_FLOAT I_1 = stress_x + stress_y + stress_z;
    MPM_FLOAT I_2 = -0.5*(I_1*I_1 - (stress_x*stress_x + stress_y*stress_y + stress_z*stress_z
        + 2.0*(_deviatoric_stress[3]*_deviatoric_stress[3] + 
               _deviatoric_stress[4]*_deviatoric_stress[4] + 
               _deviatoric_stress[5]*_deviatoric_stress[5])));
    MPM_FLOAT I_3 = stress_x*stress_y*stress_z + 
        2.0*_deviatoric_stress[3]*_deviatoric_stress[4]*_deviatoric_stress[5]
        - stress_x*_deviatoric_stress[3]*_deviatoric_stress[3] 
        - stress_y*_deviatoric_stress[4]*_deviatoric_stress[4] 
        - stress_z*_deviatoric_stress[5]*_deviatoric_stress[5];
    
    Array3D result;
    CubicFunctionRoots(1.0, -I_1, -I_2, -I_3, result);
    return move(result);
}

void PhysicalProperty::DeviatoricStressMultiplyScalar(MPM_FLOAT scalar)
{
    for (auto& element : _deviatoric_stress)
        element *= scalar;
}

void PhysicalProperty::StressRotationJaumann(SymTensor& vortex)
{
    SymTensor sigma, rotated;
    Array3D q;

    sigma[0] = _deviatoric_stress[0] + _mean_stress;
    sigma[1] = _deviatoric_stress[1] + _mean_stress;
    sigma[2] = _deviatoric_stress[2] + _mean_stress;
    sigma[3] = _deviatoric_stress[3];
    sigma[4] = _deviatoric_stress[4];
    sigma[5] = _deviatoric_stress[5];

    q[0] = 2.0*sigma[5]*vortex[2];
    q[1] = 2.0*sigma[4]*vortex[1];
    q[2] = 2.0*sigma[3]*vortex[0];

    rotated[0] = - q[0] + q[1];
    rotated[1] = + q[0] - q[2];
    rotated[2] = - q[1] + q[2];
    rotated[3] = vortex[0]*(sigma[1] - sigma[2]) + vortex[2]*sigma[4] - vortex[1]*sigma[5];
    rotated[4] = vortex[1]*(sigma[2] - sigma[0]) + vortex[0]*sigma[5] - vortex[2]*sigma[3];
    rotated[5] = vortex[2]*(sigma[0] - sigma[1]) + vortex[1]*sigma[3] - vortex[0]*sigma[4];

    for (int i = 0; i < 6; i++)
        sigma[i] += rotated[i];
    
    _mean_stress = (sigma[0] + sigma[1] + sigma[2])/3.0;

    _deviatoric_stress[0] = sigma[0] - _mean_stress;
    _deviatoric_stress[1] = sigma[1] - _mean_stress;
    _deviatoric_stress[2] = sigma[2] - _mean_stress;
    _deviatoric_stress[3] = sigma[3];
    _deviatoric_stress[4] = sigma[4];
    _deviatoric_stress[5] = sigma[5];
}

void PhysicalProperty::EquivalentStress()
{
    MPM_FLOAT J2 = 0.5*(_deviatoric_stress[0]*_deviatoric_stress[0] + 
                        _deviatoric_stress[1]*_deviatoric_stress[1] +
                        _deviatoric_stress[2]*_deviatoric_stress[2]) +
                    _deviatoric_stress[3]*_deviatoric_stress[3] +
                    _deviatoric_stress[4]*_deviatoric_stress[4] +
                    _deviatoric_stress[5]*_deviatoric_stress[5];
    _equivalent_stress = sqrt(J2*3.0);
}