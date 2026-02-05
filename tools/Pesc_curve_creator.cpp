#include <iostream>
#include <cstdint>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>

namespace cdouble
{
    double inv_2pow24 = 1.0f / 16777216.0f;
    double PI = 3.14159265f;
}

std::vector<double> loglinspace(double start, double end, int N)
{
    std::vector<double> arr_result(N);

    double log_start = log10(start);
    double log_end = log10(end);


    double log_step_size = (log_end - log_start)/(N - 1);

    for (size_t i = 0; i < N; i++)
    {
        arr_result[i] = pow(10, log_start + (static_cast<double>(i) * log_step_size));
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

    inline double next_double() {
        // Take upper 24 bits of next() and convert to double in [0,1)
        return (next() >> 40) * cdouble::inv_2pow24;
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

struct Vec3double
{
    double x, y, z;
    Vec3double(double x_in, double y_in, double z_in) : x(x_in), y(y_in), z(z_in) {}
};


Vec3double generate_angle_HG_double(double dx, double dy, double dz, double g, FastRNG& rng)
{
    double mu;
    if (fabsf(g) < 1e-6f) {
        mu = rng.uniform()*2.f - 1.f;
    } else {
        double xi = rng.uniform();
        double gg = g * g;
        double temp = (1.f - gg) / (1.f - g + 2.f * g * xi);
        mu = (1.f + gg - temp * temp) / (2.f * g);
    }

    double phi = rng.uniform()*2*cdouble::PI;

    double sin_phi = sinf(phi);
    double cos_phi = cosf(phi);
    double sin_theta = sqrtf(1.f - mu*mu);

    // Creating tangent and bi-tangent lines (t and b)
    double tdx, tdy, tdz, bdx, bdy, bdz, vdx, vdy, vdz;
    double dx2dy2 = dx*dx + dy*dy;
    double inv_tnorm = 1.f/sqrtf(dx2dy2);
    double inv_bnorm = 1.f/sqrtf(dx2dy2 + dz*dz);

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

    Vec3double vec_out(vdx, vdy, vdz);

    return vec_out;
}


double f_Pesc_MC3D(const double cell_dx,
                  const double cell_dy,
                  const double cell_dz,
                  const double cell_kext,
                  const int N)
{
    const int Nbatch = 1024;
    const double x0 = 0.f;
    const double y0 = 0.f;
    const double z0 = 0.f;
    const double x1 = x0 + cell_dx;
    const double y1 = y0 + cell_dy;
    const double z1 = z0 + cell_dz;

    int Nout = 0;

    const double tiny_eps = 1e-12;
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        FastRNG rng_local(std::chrono::high_resolution_clock::now().time_since_epoch().count() + tid);
        int Nout_local = 0;

        #pragma omp for schedule(dynamic, Nbatch)
        for (long int idx_photon = 0; idx_photon < N; idx_photon++)
        {
            // Generating the photon
            double x = rng_local.uniform()*cell_dx;
            double y = rng_local.uniform()*cell_dy;
            double z = rng_local.uniform()*cell_dz;

            double mu = rng_local.uniform()*2.f - 1.f;
            double az = rng_local.uniform()*2.f * 2*cdouble::PI;

            double s = sqrtf(1.f - mu*mu);
            double dx = s*cosf(az);
            double dy = s*sinf(az);
            double dz = mu;
            
            dx = copysignf(fmaxf(fabsf(dx), tiny_eps), dx);
            dy = copysignf(fmaxf(fabsf(dy), tiny_eps), dy);
            dz = copysignf(fmaxf(fabsf(dz), tiny_eps), dz);

            double tau = -logf(1.f - rng_local.uniform());

            // Propagating
            double f = tau / cell_kext;
            double x_new = x + dx*f;
            double y_new = y + dy*f;
            double z_new = z + dz*f;

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

    double Pesc = ((double) Nout) / N;
    return Pesc;
}



void write_Pesc_curve(double dx, double dy, double dz, int Nphot, int N_curvepoints)
{
    double kext_low = 1e-10;
    double kext_high = 1e4;
    
    std::vector<double> arr_kext = loglinspace(kext_low, kext_high, N_curvepoints);

    std::ostringstream filename;
    filename << "../data_input/Esc_curves/Esc_kext_curve_"
             << std::fixed << std::setprecision(0) << dx << "_"
             << std::fixed << std::setprecision(0) << dy << "_"
             << std::fixed << std::setprecision(2) << dz << ".csv";
    std::ofstream file(filename.str());
    if (!file.is_open())
    {
        std::cout << "Could not open " + filename.str() << "!" << std::endl;
    }
    // Writing Header
    file << "kext,Esc\n";

    file << std::scientific << std::setprecision(15);

    for (int i = 0; i < N_curvepoints; i++)
    {
        double cell_kext = arr_kext[i];
        double Pesc = f_Pesc_MC3D(dx, dy, dz, cell_kext, Nphot);
        file << cell_kext << ',' << Pesc << "\n";
        std::cout << dx << ',' << dy << ',' << dz << ',' << i << std::endl;
    }
}


int main(int argc, char* argv[])
{
    // Arg 1: dx
    // Arg 2: dy
    // Arg 3: dz
    // Arg 4: nphot_pow
    // Arg 5: N_curvepoints
    
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;



    int nphot_pow = std::stoi(argv[4]);


    int Nphot = pow(2, nphot_pow);

    double dx = std::stod(argv[1]);
    double dy = std::stod(argv[2]);
    double dz = std::stod(argv[3]);
    int N_curvepoints = std::stoi(argv[5]);
    

    std::cout << "Starting..." << std::endl;
    write_Pesc_curve(dx, dy, dz, Nphot, N_curvepoints);

    return 0;
}