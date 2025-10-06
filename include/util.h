#pragma once

#include <vector>
#include <iostream>

template <typename T>
void LOGvec(const std::vector<T>& arr, const std::string& name = "Array", bool enters = false)
{
    std::cout << name << ": [";
    for (size_t i = 0; i < arr.size(); i++) 
    {
        std::cout << arr[i];
        
        if (i < (arr.size()-1))
        {
            std::cout << ", ";
        } else
        {
            std::cout << ']' << std::endl;
        }
        
        if (enters)
        {
            std::cout << std::endl;
        }
    }
}


namespace cf
{
    const float PI = 3.1415927;
    const float RHO = 1.225;
    const float CP = 1004;
}
