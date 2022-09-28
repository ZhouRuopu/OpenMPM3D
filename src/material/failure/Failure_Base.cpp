#include "./Failure_Base.h"
Failure_Base::Failure_Base()
{
    Type = "";
    Erosion = false;
}

Failure_Base::~Failure_Base()
{
}

bool Failure_Base::Initialize(map<string, MPM_FLOAT> &failure_para)
{
    for(map<string, MPM_FLOAT>::iterator iter = failure_para.begin(); 
        iter != failure_para.end(); iter++)
    {
        if(ParameterMap_Failure.find(iter->first) != ParameterMap_Failure.end())
            *ParameterMap_Failure[iter->first] = iter->second;
        else
        {
            string error_msg = "Can't find the failure model parameter " + iter->first + " at " + Type;
            MPM3D_ErrorMessage(error_msg);
            return false;
        }
    }
    return true;
}