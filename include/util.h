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

inline double kahan_sum(const std::vector<double>& values) 
{
    // from Mr. Chat
    double sum = 0.0;
    double c = 0.0;  // Compensation
    for (double x : values) {
        double y = x - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}


namespace cf
{
    const double PI = 3.1415927;
    const double RHO = 1.225;
    const double CP = 1004;
}
