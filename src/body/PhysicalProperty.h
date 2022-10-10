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
    Info: Class definition for physical properties
    Code-writter: Ruichen Ni
    Date: 2022.9.29
==============================================================*/

#ifndef _PhysicalProperty_H_
#define _PhysicalProperty_H_
#include "../main/MPM3D_MACRO.h"

class PhysicalProperty
{
public:
    PhysicalProperty();
    PhysicalProperty(const PhysicalProperty& pp);
    ~PhysicalProperty();

    //!> Allocate Memory For Extra Particle Property
    void AllocateMemoryForExtraParticleProperty(int number);

    //!> override operator [] to get extra particle property
    inline MPM_FLOAT& operator[] (int index)
    {
        return _extra_properties[_extra_property_positions[index]];
    }

    //!> Update the volume based on the incremental volumetric strain
    inline void UpdateVolume(MPM_FLOAT (&de)[6])
    {
        _volume *= (1 + de[0] + de[1] + de[2]);
        _density = _mass/_volume;
    }

    inline void UpdateVolume_Exponent(MPM_FLOAT (&de)[6])
    {
        _volume *= exp(de[0] + de[1] + de[2]);
        _density = _mass/_volume;
    }

    //!> Calculate the principle stress
    Array3D&& CalculatePrincipleStress();

    //!> Multiply the deviatoric stress with a scalar
    void DeviatoricStressMultiplyScalar(MPM_FLOAT scalar);

    //!> Rotate the stress with Jaumann rate
    void StressRotationJaumann(SymTensor& vortex);

    //!> Calculate Equivalent Stress
    void EquivalentStress();
private:
    MPM_FLOAT _mass;
    MPM_FLOAT _volume;
    MPM_FLOAT _density;
    MPM_FLOAT _mean_stress;
    SymTensor _deviatoric_stress;   //!< SDxx, SDyy, SDzz, SDyz, SDxz, SDxy in sequence
    MPM_FLOAT _equivalent_stress;   //!< Von Mises Stress
    MPM_FLOAT _bulk_q;              //!< artificial bulk viscosity
    MPM_FLOAT _internal_energy;
    MPM_FLOAT _sound_speed;
    bool _failure;
    bool _eroded;

    //!> Extra Particle Properties
    MPM_FLOAT* _extra_properties;
    int* _extra_property_positions;

public:
//!> various Get/Set function
    inline MPM_FLOAT GetMass() {return _mass;};
    inline void SetMass(MPM_FLOAT mass) {_mass = mass;};

    inline MPM_FLOAT GetVolume() {return _volume;};
    inline void SetVolume(MPM_FLOAT vol) {_volume = vol;};

    inline MPM_FLOAT GetDensity() {return _density;};
    inline void SetDensity(MPM_FLOAT rho) {_density = rho;};

    inline MPM_FLOAT GetMeanStress() {return _mean_stress;};
    inline void SetMeanStress(MPM_FLOAT SM) {_mean_stress = SM;};

    inline SymTensor GetDeviatoricStress() {return _deviatoric_stress;}
    inline void SetDeviatoricStress(SymTensor& sd) {_deviatoric_stress = sd;}

    inline MPM_FLOAT GetEquivalentStress() {return _equivalent_stress;}
    inline void SetEquivalentStress(MPM_FLOAT seqv) {_equivalent_stress = seqv;}

    inline MPM_FLOAT GetBulkViscosity() {return _bulk_q;};
    inline void SetBulkViscosity(MPM_FLOAT q) {_bulk_q = q;};

    inline MPM_FLOAT GetInternalEnergy() {return _internal_energy;};
    inline void SetInternalEnergy(MPM_FLOAT ie) {_internal_energy = ie;};

    inline MPM_FLOAT GetSoundSpeed() {return _sound_speed;};
    inline void SetSoundSpeed(MPM_FLOAT c) {_sound_speed = c;};

    inline bool is_Failed() {return _failure;};
    inline void Failed() {_failure = true;};

    inline bool is_Eroded() {return _eroded;};
    inline void Eroded() {_eroded = true;};
};

#endif