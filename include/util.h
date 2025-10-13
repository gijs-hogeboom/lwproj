#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>

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

template<typename T>
void printElement(const T& value, int width = 10) {
    std::cout << std::setw(width) << value;
}

template<typename FirstVec, typename... RestVecs>
void LOGvecCompare(const FirstVec& first, const RestVecs&... rest) {
    // Ensure all vectors have the same size
    size_t n = first.size();
    ((rest.size() == n ? void() : throw std::runtime_error("All vectors must have the same size")), ...);

    // Print each row
    for (size_t i = 0; i < n; ++i) {
        printElement(first[i]);
        ((std::cout << ", ", printElement(rest[i])), ...);
        std::cout << "," << std::endl;
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
