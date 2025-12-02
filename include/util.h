#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cstdint>
#include <queue>
#include <numeric>
#include <algorithm>
#include <fstream>


// Constants
namespace cfloat
{
    const float PI = 3.1415927;
    const float RHO = 1.225;
    const float CP = 1004;
    const float u = pow(2, -24);
}

namespace cdouble
{
    const double PI = 3.141592653589793;
    const double RHO = 1.225;
    const double CP = 1004.;
    const double u = pow(2, -53);
}




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


inline float estimate_max_numerical_float_error(float s)
{
    return cfloat::u * s;
}


inline double estimate_max_numerical_double_error(double s)
{
    return cdouble::u * s;
}


inline std::vector<double> linspace(double start, double end, int N)
{
    std::vector<double> arr_result(N);

    double step_size = (end - start)/(N - 1);

    for (size_t i = 0; i < N; i++)
    {
        arr_result[i] = start + (static_cast<double>(i) * step_size);
    }

    return arr_result;
}




inline double trapezoid(const std::vector<double>& arr_x,
                const std::vector<double>& arr_y)
{
    size_t n = arr_x.size();
    std::vector<double> values((n - 1));

    for (size_t i = 0; i < (n - 1); i++)
    {
        double dx = arr_x[i+1] - arr_x[i];
        double Y  = arr_y[i] * dx + (arr_y[i+1] - arr_y[i]) * dx / 2.0;
        values[i] = Y;
    }

    double result = std::accumulate(values.begin(), values.end(), 0.0);
    return result;
}





struct Xoshiro256ss {
    uint64_t s[4];

    Xoshiro256ss(uint64_t seed = 1) {
        // SplitMix64 to initialize the state
        uint64_t z = seed + 0x9e3779b97f4a7c15ULL;
        for (int i = 0; i < 4; ++i) {
            z ^= (z >> 30); z *= 0xbf58476d1ce4e5b9ULL;
            z ^= (z >> 27); z *= 0x94d049bb133111ebULL;
            z ^= (z >> 31);
            s[i] = z;
            z += 0x9e3779b97f4a7c15ULL;
        }
    }

    inline uint64_t next() {
        const uint64_t result = rotl(s[1] * 5, 7) * 9;
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }

    inline double next_double() {
        // Take upper 53 bits of next() and convert to double in [0,1)
        return (next() >> 11) * (1.0 / 9007199254740992.0);
    }

    inline float next_float() {
        // Take upper 24 bits of next() and convert to float in [0,1)
        return (next() >> 40) * (1.0 / 16777216.0);
    }

private:
    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
};


struct FastRNG {
    Xoshiro256ss rng;
    FastRNG(uint64_t seed) : rng(seed) {}

    inline double uniform() { return rng.next_double(); }

    inline int uniform_int(int max_exclusive) {
        return static_cast<int>(((unsigned __int128)rng.next() * (unsigned __int128)max_exclusive) >> 64);
    }
};


struct AliasTable_double {
    std::vector<double> prob;
    std::vector<int> alias;
    std::vector<double> weights;
    int n;

    AliasTable_double(const std::vector<double>& weights_in) {
        n = weights_in.size();
        prob.resize(n);
        alias.resize(n);

        std::vector<double> scaled(weights_in);
        double sum = std::accumulate(scaled.begin(), scaled.end(), 0.0);
        for (auto& w : scaled) w *= n / sum;

        weights.resize(n);
        for (int i = 0; i < n; i++) weights[i] = weights_in[i] / sum;
        


        std::queue<int> small, large;
        for (int i = 0; i < n; ++i)
            (scaled[i] < 1.0 ? small : large).push(i);

        while (!small.empty() && !large.empty()) {
            int s = small.front(); small.pop();
            int l = large.front(); large.pop();
            prob[s] = scaled[s];
            alias[s] = l;
            scaled[l] = scaled[l] + scaled[s] - 1.0;
            (scaled[l] < 1.0 ? small : large).push(l);
        }

        while (!large.empty()) { prob[large.front()] = 1.0; large.pop(); }
        while (!small.empty()) { prob[small.front()] = 1.0; small.pop(); }
    }

    inline int sample(FastRNG& rng) const {
        int i = rng.rng.next() % n;
        double r = rng.rng.next_double();
        return (r < prob[i]) ? i : alias[i];
    }
};


struct AliasTable_float {
    std::vector<float> prob;
    std::vector<int> alias;
    std::vector<float> weights;
    int n;

    AliasTable_float(const std::vector<float>& weights_in) {
        n = weights_in.size();
        prob.resize(n);
        alias.resize(n);

        std::vector<float> scaled(weights_in);
        float sum = std::accumulate(scaled.begin(), scaled.end(), 0.0);
        for (auto& w : scaled) w *= n / sum;

        weights.resize(n);
        for (int i = 0; i < n; i++) weights[i] = weights_in[i] / sum;

        std::queue<int> small, large;
        for (int i = 0; i < n; ++i)
            (scaled[i] < 1.0 ? small : large).push(i);

        while (!small.empty() && !large.empty()) {
            int s = small.front(); small.pop();
            int l = large.front(); large.pop();
            prob[s] = scaled[s];
            alias[s] = l;
            scaled[l] = scaled[l] + scaled[s] - 1.0;
            (scaled[l] < 1.0 ? small : large).push(l);
        }

        while (!large.empty()) { prob[large.front()] = 1.0; large.pop(); }
        while (!small.empty()) { prob[small.front()] = 1.0; small.pop(); }
    }

    inline int sample(FastRNG& rng) const {
        int i = rng.rng.next() % n;
        float r = rng.rng.next_float();
        return (r < prob[i]) ? i : alias[i];
    }
};



inline std::string f_dz_string(double value) {
    std::ostringstream temp;
    temp << std::fixed << std::setprecision(5) << value;
    std::string out = temp.str();
    out = out.substr(0, out.size() - 3);
    return out;
}

inline std::string f_Pesccurve_name(int dx, int dy, double dz)
{
    std::string dzs = f_dz_string(dz);
    std::ostringstream out;
    out << "/home/gijs-hogeboom/dev/lwproj/data_input/Esc_curves/Esc_kext_curve_" << (int) dx << "_" << (int) dy << "_" << dzs << ".csv";
    std::string filename = out.str();
    return filename;
}






class LinearInterpolator_double {
public:
    LinearInterpolator_double(const std::vector<double>& x,
                              const std::vector<double>& y)
        : xs(x), ys(y)
    {
        if (xs.size() != ys.size() || xs.size() < 2) {
            throw std::invalid_argument("x and y arrays must have same size >= 2");
        }

        // ensure strictly increasing x
        for (size_t i = 1; i < xs.size(); ++i) {
            if (xs[i] <= xs[i - 1]) {
                throw std::invalid_argument("x values must be strictly increasing");
            }
        }
    }

    // interpolate y at value xq
    inline double operator()(double xq) const {
        // out-of-range -> throw
        if (xq < xs.front() || xq > xs.back()) {
            throw std::out_of_range("query x is outside interpolation range");
        }

        // find first element greater than xq
        auto it = std::lower_bound(xs.begin(), xs.end(), xq);
        
        if (it == xs.begin())
            return ys.front();
        
        if (it == xs.end())
            return ys.back();

        // indices of bounding interval
        size_t i1 = it - xs.begin();
        size_t i0 = i1 - 1;

        double x0 = xs[i0], x1 = xs[i1];
        double y0 = ys[i0], y1 = ys[i1];

        // linear interpolation
        double t = (xq - x0) / (x1 - x0);
        return y0 + t * (y1 - y0);
    }

private:
    std::vector<double> xs;
    std::vector<double> ys;
};



inline std::size_t count_lines(std::fstream& file) {
    std::size_t count = 0;
    std::string line;

    while (std::getline(file, line))
        ++count;

    return count;
}
