#include <iostream>
#include <cstdint>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>

namespace cfloat
{
    float inv_2pow24 = 1.0f / 16777216.0f;
    float PI = 3.14159265f;
}

std::vector<float> loglinspace(float start, float end, int N)
{
    std::vector<float> arr_result(N);

    float log_start = log10(start);
    float log_end = log10(end);


    float log_step_size = (log_end - log_start)/(N - 1);

    for (size_t i = 0; i < N; i++)
    {
        arr_result[i] = pow(10, log_start + (static_cast<float>(i) * log_step_size));
    }

    return arr_result;
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

    // inline double next_double() {
    //     // Take upper 53 bits of next() and convert to double in [0,1)
    //     return (next() >> 11) * (1.0 / 9007199254740992.0);
    // }

    inline float next_float() {
        // Take upper 24 bits of next() and convert to float in [0,1)
        return (next() >> 40) * cfloat::inv_2pow24;
    }

private:
    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
};


struct FastRNG {
    Xoshiro256ss rng;
    FastRNG(uint64_t seed) : rng(seed) {}

    inline float uniform() { return rng.next_float(); }

    inline int uniform_int(int max_exclusive) {
        return static_cast<int>(((unsigned __int128)rng.next() * (unsigned __int128)max_exclusive) >> 64);
    }
};

struct Vec3float
{
    float x, y, z;
    Vec3float(float x_in, float y_in, float z_in) : x(x_in), y(y_in), z(z_in) {}
};


Vec3float generate_angle_HG_float(float dx, float dy, float dz, float g, FastRNG& rng)
{
    float mu;
    if (fabsf(g) < 1e-6f) {
        mu = rng.uniform()*2.f - 1.f;
    } else {
        float xi = rng.uniform();
        float gg = g * g;
        float temp = (1.f - gg) / (1.f - g + 2.f * g * xi);
        mu = (1.f + gg - temp * temp) / (2.f * g);
    }

    float phi = rng.uniform()*2*cfloat::PI;

    float sin_phi = sinf(phi);
    float cos_phi = cosf(phi);
    float sin_theta = sqrtf(1.f - mu*mu);

    // Creating tangent and bi-tangent lines (t and b)
    float tdx, tdy, tdz, bdx, bdy, bdz, vdx, vdy, vdz;
    float dx2dy2 = dx*dx + dy*dy;
    float inv_tnorm = 1.f/sqrtf(dx2dy2);
    float inv_bnorm = 1.f/sqrtf(dx2dy2 + dz*dz);

    // t = up "x" u, normalized to length = 1
    // b = u "x" t, normalized to length = 1

    tdx = -dy * inv_tnorm;
    tdy = dx * inv_tnorm;
    tdz = 0.f;

    bdx = -tdy * dz * inv_bnorm;
    bdy = tdx * dz * inv_bnorm;
    bdz = dx2dy2 * inv_tnorm * inv_bnorm;

    vdx =   sin_theta * cos_phi * tdx +   sin_theta * sin_phi * bdx + mu * dx;
    vdy =   sin_theta * cos_phi * tdy +   sin_theta * sin_phi * bdy + mu * dy;
    vdz = /*sin_theta * cos_phi * tdz +*/ sin_theta * sin_phi * bdz + mu * dz; // tdz = 0.

    Vec3float vec_out(vdx, vdy, vdz);

    return vec_out;
}


float f_Pesc_MC3D(const float cell_dx,
                  const float cell_dy,
                  const float cell_dz,
                  const float cell_kext,
                  const int N)
{
    const int Nbatch = 1024;
    const float x0 = 0.f;
    const float y0 = 0.f;
    const float z0 = 0.f;
    const float x1 = x0 + cell_dx;
    const float y1 = y0 + cell_dy;
    const float z1 = z0 + cell_dz;

    int Nout = 0;

    const float tiny_eps = 1e-12;
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        FastRNG rng_local(std::chrono::high_resolution_clock::now().time_since_epoch().count() + tid);
        int Nout_local = 0;

        #pragma omp for schedule(dynamic, Nbatch)
        for (long int idx_photon = 0; idx_photon < N; idx_photon++)
        {
            // Generating the photon
            float x = rng_local.uniform()*cell_dx;
            float y = rng_local.uniform()*cell_dy;
            float z = rng_local.uniform()*cell_dz;

            float mu = rng_local.uniform()*2.f - 1.f;
            float az = rng_local.uniform()*2.f * 2*cfloat::PI;

            float s = sqrtf(1.f - mu*mu);
            float dx = s*cosf(az);
            float dy = s*sinf(az);
            float dz = mu;
            
            dx = copysignf(fmaxf(fabsf(dx), tiny_eps), dx);
            dy = copysignf(fmaxf(fabsf(dy), tiny_eps), dy);
            dz = copysignf(fmaxf(fabsf(dz), tiny_eps), dz);

            float tau = -logf(1.f - rng_local.uniform());

            // Propagating
            float f = tau / cell_kext;
            float x_new = x + dx*f;
            float y_new = y + dy*f;
            float z_new = z + dz*f;

            // Tallying if photon escaped
            if ((x_new < x0) || (x_new > x1) ||
                (y_new < y0) || (y_new > y1) ||
                (z_new < z0) || (z_new > z1))
            {
                Nout_local += 1;   
            }
        }

        #pragma omp critical
        Nout += Nout_local;
    }

    float Pesc = ((float) Nout) / N;
    return Pesc;
}



void write_Pesc_curve(float dx, float dy, float dz, int Nphot)
{
    float kext_low = 1e-10;
    float kext_high = 1e4;
    int N_curvepoints = 400;

    std::vector<float> arr_kext = loglinspace(kext_low, kext_high, N_curvepoints);

    std::ostringstream filename;
    filename << "../data_input/Esc_curves_new/Esc_kext_curve_"
             << std::fixed << std::setprecision(2) << dx << "_"
             << std::fixed << std::setprecision(2) << dy << "_"
             << std::fixed << std::setprecision(2) << dz << ".csv";
    std::ofstream file(filename.str());
    if (!file.is_open())
    {
        std::cout << "Could not open " + filename.str() << "!" << std::endl;
    }
    // Writing Header
    file << "kext,Esc\n";

    for (int i = 0; i < N_curvepoints; i++)
    {
        float cell_kext = arr_kext[i];
        float Pesc = f_Pesc_MC3D(dx, dy, dz, cell_kext, Nphot);
        file << cell_kext << ',' << Pesc << "\n";
        std::cout << i << std::endl;
    }
}


int main(int argc, char* argv[])
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;



    int nphot_pow = 25;


    int Nphot = pow(2, nphot_pow);

    float dx = 1.f;
    float dy = 1.f;
    float dz = 1.f;
    float kext = 1.f;

    std::cout << "Starting..." << std::endl;
    write_Pesc_curve(dx, dy, dz, Nphot);

    return 0;
}