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

inline float kahan_sum(const std::vector<float>& values) 
{
    // from Mr. Chat
    float sum = 0.0;
    float c = 0.0;  // Compensation
    for (float x : values) {
        float y = x - c;
        float t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}


namespace cf
{
    const float PI = 3.1415927;
    const float RHO = 1.225;
    const float CP = 1004;
}
