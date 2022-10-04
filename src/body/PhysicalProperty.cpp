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

PhysicalProperty::PhysicalProperty()
{
    _mass = 0.0;
    _volume = 0.0;
    _mean_stress = 0.0;
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
    _mean_stress = pp._mean_stress;
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