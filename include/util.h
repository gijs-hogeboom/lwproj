#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cstdint>
#include <cmath>
#include <queue>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <limits>

#include "precision.h"



// Constants
namespace constants
{
    const Real PI = 3.14159265358979323846264338;
    const Real RHO = 1.225;
    const Real CP = 1004;
    const Real u = pow(2, -24);
    const Real inv_2pow24 = 1.0f / 16777216.0f;
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

template<typename FirstVal, typename... RestVals>
void LOGvars(const FirstVal first, const RestVals... rest)
{
    printElement(first);
    ((std::cout << " | ", printElement(rest)), ...);
    std::cout << std::endl;
}


inline Real kahan_sum(const std::vector<Real>& values) 
{
    // from Mr. Chat
    Real sum = 0.0;
    Real c = 0.0;  // Compensation
    for (Real x : values) {
        Real y = x - c;
        Real t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}


inline std::vector<Real> linspace(Real start, Real end, int N)
{
    std::vector<Real> arr_result(N);

    Real step_size = (end - start)/(N - 1);

    for (size_t i = 0; i < N; i++)
    {
        arr_result[i] = start + (static_cast<Real>(i) * step_size);
    }

    return arr_result;
}



inline Real trapezoid(const std::vector<Real>& arr_x,
                             const std::vector<Real>& arr_y)
{
    size_t n = arr_x.size();
    std::vector<Real> values((n - 1));

    for (size_t i = 0; i < (n - 1); i++)
    {
        Real dx = arr_x[i+1] - arr_x[i];
        Real Y  = arr_y[i] * dx + (arr_y[i+1] - arr_y[i]) * dx / static_cast<Real>(2.);
        values[i] = Y;
    }

    Real result = std::accumulate(values.begin(), values.end(), static_cast<Real>(0.));
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

    inline Real next_real()
    {
        constexpr int mantissa_bits = std::numeric_limits<Real>::digits;
        constexpr Real inv = Real(1.) / Real(std::uint64_t(1) << mantissa_bits);

        return Real(next() >> (64 - mantissa_bits)) * inv;
    }

private:
    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
};


struct FastRNG {
    Xoshiro256ss rng;
    FastRNG(uint64_t seed) : rng(seed) {}

    inline Real uniform()
    { 
        return rng.next_real(); 
    }

    inline int uniform_int(int max_exclusive) 
    {
        return static_cast<int>(((unsigned __int128)rng.next() * (unsigned __int128)max_exclusive) >> 64);
    }
};



struct AliasTable {
    std::vector<Real> prob;
    std::vector<int> alias;
    std::vector<Real> weights;
    int n;

    AliasTable(const std::vector<Real>& weights_in) {
        n = weights_in.size();
        prob.resize(n);
        alias.resize(n);

        std::vector<Real> scaled(weights_in);
        Real sum = std::accumulate(scaled.begin(), scaled.end(), static_cast<Real>(0.));
        for (auto& w : scaled) w *= n / sum;

        weights.resize(n);
        for (int i = 0; i < n; i++) weights[i] = weights_in[i] / sum;

        std::queue<int> small, large;
        for (int i = 0; i < n; ++i)
            (scaled[i] < static_cast<Real>(1.) ? small : large).push(i);

        while (!small.empty() && !large.empty()) {
            int s = small.front(); small.pop();
            int l = large.front(); large.pop();
            prob[s] = scaled[s];
            alias[s] = l;
            scaled[l] = scaled[l] + scaled[s] - static_cast<Real>(1.);
            (scaled[l] < static_cast<Real>(1.) ? small : large).push(l);
        }

        while (!large.empty()) { prob[large.front()] = 1.0; large.pop(); }
        while (!small.empty()) { prob[small.front()] = 1.0; small.pop(); }
    }

    inline int sample(FastRNG& rng) const {
        int i = rng.rng.next() % n;
        float r = rng.rng.next_real();
        return (r < prob[i]) ? i : alias[i];
    }
};


inline std::size_t count_lines(std::fstream& file) {
    std::size_t count = 0;
    std::string line;

    while (std::getline(file, line))
        ++count;

    return count;
}

inline std::string f_dz_string(Real value) {
    std::ostringstream temp;
    temp << std::fixed << std::setprecision(2) << value;
    std::string out = temp.str();
    return out;
}

inline std::string f_Pesccurve_name(int dx, int dy, Real dz)
{
    std::string dzs = f_dz_string(dz);
    std::ostringstream out;
    out << "../data_input/Esc_curves/Esc_kext_curve_" << (int) dx << "_" << (int) dy << "_" << dzs << ".csv";
    std::string filename = out.str();
    return filename;
}



class LinearInterpolator {
public:
    LinearInterpolator(const std::vector<Real>& x,
                       const std::vector<Real>& y)
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
    inline Real operator()(Real xq) const {
        // out-of-range -> throw
        if (xq < xs.front() || xq > xs.back()) {
            throw std::out_of_range("query x (" + std::to_string(xq) + ") is outside interpolation range");
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

        Real x0 = xs[i0], x1 = xs[i1];
        Real y0 = ys[i0], y1 = ys[i1];

        // linear interpolation
        Real t = (xq - x0) / (x1 - x0);
        return y0 + t * (y1 - y0);
    }

private:
    std::vector<Real> xs;
    std::vector<Real> ys;
};




struct Vec3
{
    Real x, y, z;
    Vec3(Real x_in, Real y_in, Real z_in) : x(x_in), y(y_in), z(z_in) {}
};

inline Vec3 generate_angle_HG(Real dx, Real dy, Real dz, Real g, FastRNG& rng)
{
    Real mu;
    if (fabsf(g) < static_cast<Real>(1e-6)) {
        mu = rng.uniform()*static_cast<Real>(2.) - static_cast<Real>(1.);
    } else {
        Real xi = rng.uniform();
        Real gg = g * g;
        Real temp = (static_cast<Real>(1.) - gg) / (static_cast<Real>(1.) - g + static_cast<Real>(2.) * g * xi);
        mu = (static_cast<Real>(1.) + gg - temp * temp) / (static_cast<Real>(2.) * g);
    }

    Real phi = rng.uniform()*static_cast<Real>(2.)*constants::PI;

    Real sin_phi = std::sin(phi);
    Real cos_phi = std::cos(phi);
    Real sin_theta = std::sqrt(static_cast<Real>(1.) - mu*mu);

    // Creating tangent and bi-tangent lines (t and b)
    Real tdx, tdy, tdz, bdx, bdy, bdz, vdx, vdy, vdz;
    Real dx2dy2 = dx*dx + dy*dy;
    Real inv_tnorm = static_cast<Real>(1.)/std::sqrt(dx2dy2);
    Real inv_bnorm = static_cast<Real>(1.)/std::sqrt(dx2dy2 + dz*dz);

    // t = up "x" u, normalized to length = 1
    // b = u "x" t, normalized to length = 1

    tdx = -dy * inv_tnorm;
    tdy = dx * inv_tnorm;
    tdz = 0.;

    bdx = -tdy * dz * inv_bnorm;
    bdy = tdx * dz * inv_bnorm;
    bdz = dx2dy2 * inv_tnorm * inv_bnorm;

    vdx =   sin_theta * cos_phi * tdx +   sin_theta * sin_phi * bdx + mu * dx;
    vdy =   sin_theta * cos_phi * tdy +   sin_theta * sin_phi * bdy + mu * dy;
    vdz = /*sin_theta * cos_phi * tdz +*/ sin_theta * sin_phi * bdz + mu * dz; // tdz = 0.

    Vec3 vec_out(vdx, vdy, vdz);

    return vec_out;
}   
