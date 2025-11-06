// ANCIENNE VERSION SIMPLIFIEE : uniquement value-noise !!!!!

// #pragma once
// #include <cmath>
// #include <cstdint>
// #include <algorithm>

// struct HeightParams {
//     int   N = 257;        // sommets par côté (puissance de 2 + 1)
//     float scale = 1.0f;   // pas entre points en X/Z (monde)
//     float amp = 10.0f;    // amplitude globale
//     float freq = 0.8f;    // fréquence de base
//     int   octaves = 4;    // 3..6 ok
//     float lacunarity = 2.0f;
//     float gain = 0.5f;
//     int   seed = 42;
//     bool  center = true;  // centrer le terrain autour de (0,0,0)
// };

// // --- helpers ---
// static inline float lerp_f(float a,float b,float t){ return a+(b-a)*t; }
// static inline float fade_f(float t){ return t*t*t*(t*(t*6.0f-15.0f)+10.0f); } // Perlin's fade

// // petit hash déterministe → [0,1]
// static inline float hash2i(int x, int y, int seed) {
//     uint32_t h = uint32_t(x)*374761393u + uint32_t(y)*668265263u + uint32_t(seed)*2246822519u;
//     h = (h ^ (h >> 13)) * 1274126177u;
//     return (h ^ (h >> 16)) * (1.0f/4294967295.0f);
// }

// // value-noise bilinéaire lissé, retourne [-1,1]
// static inline float valueNoise2D(float x, float y, int seed){
//     int xi = int(std::floor(x));
//     int yi = int(std::floor(y));
//     float tx = x - xi, ty = y - yi;

//     float a = hash2i(xi,   yi,   seed);
//     float b = hash2i(xi+1, yi,   seed);
//     float c = hash2i(xi,   yi+1, seed);
//     float d = hash2i(xi+1, yi+1, seed);

//     float ux = fade_f(tx), uy = fade_f(ty);
//     float v = lerp_f( lerp_f(a,b,ux), lerp_f(c,d,ux), uy );
//     return v*2.0f - 1.0f;
// }

// // fBm (Fractional Brownian Motion)
// static inline float fbm(float x, float y, const HeightParams& P){
//     float f = P.freq, a = P.amp, sum = 0.0f, norm = 0.0f;
//     for(int o=0;o<P.octaves;++o){
//         sum  += a * valueNoise2D(x*f, y*f, P.seed + 131*o);
//         norm += a;
//         f *= P.lacunarity;
//         a *= P.gain;
//     }
//     return (norm>1e-6f) ? (sum/norm) : 0.0f;
// }

// =================== Version complète avec Perlin, ridged, warping, bassins ===================

#pragma once
#include <cmath>
#include <cstdint>
#include <algorithm>

enum class NoiseType { Value, Perlin, PerlinRidged };

struct HeightParams {
    int   N = 257;
    float scale = 1.0f;
    float amp = 10.0f;
    float freq = 0.8f;
    int   octaves = 4;
    float lacunarity = 2.0f;
    float gain = 0.5f;
    int   seed = 42;
    bool  center = true;

    NoiseType model = NoiseType::Perlin;
    bool  useWarp   = true;
    bool  useBasins = true;
    float ridgedPow = 2.0f;
    float warpStrength = 0.0f;
    float warpFreq = 0.0f;
    float heightBias = 0.0f;

    HeightParams() = default;
    HeightParams(int N_, float scale_, float amp_, float freq_, int oct_, float lac, float g,
                 int seed_, bool center_, bool ridged_=false, float warpS_=0.f, float warpF_=0.f)
        : N(N_), scale(scale_), amp(amp_), freq(freq_), octaves(oct_), lacunarity(lac),
          gain(g), seed(seed_), center(center_),
          model(ridged_ ? NoiseType::PerlinRidged : NoiseType::Perlin),
          useWarp(warpS_>0.f), useBasins(false), ridgedPow(2.0f),
          warpStrength(warpS_), warpFreq(warpF_), heightBias(0.0f) {}
};



// =================== Helpers ===================
static inline float lerp_f(float a, float b, float t) { return a + (b - a) * t; }
static inline float fade_f(float t) { return t*t*t*(t*(t*6.0f - 15.0f) + 10.0f); }

static inline uint32_t hash2_32(int x, int y, int seed) {
    uint32_t h = uint32_t(x) * 374761393u
               + uint32_t(y) * 668265263u
               + uint32_t(seed) * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

static inline void grad2_from_hash(uint32_t h, float& gx, float& gy) {
    static const float G[8][2] = {
        { 1.0f, 0.0f}, { 0.70710678f,  0.70710678f},
        { 0.0f, 1.0f}, {-0.70710678f,  0.70710678f},
        {-1.0f, 0.0f}, {-0.70710678f, -0.70710678f},
        { 0.0f,-1.0f}, { 0.70710678f, -0.70710678f}
    };
    const int i = int(h & 7u);
    gx = G[i][0]; gy = G[i][1];
}

// =================== Perlin 2D ===================
static inline float perlin2D(float x, float y, int seed) {
    int xi = int(std::floor(x));
    int yi = int(std::floor(y));
    float tx = x - xi;
    float ty = y - yi;

    float g00x, g00y, g10x, g10y, g01x, g01y, g11x, g11y;
    grad2_from_hash(hash2_32(xi,   yi,   seed), g00x, g00y);
    grad2_from_hash(hash2_32(xi+1, yi,   seed), g10x, g10y);
    grad2_from_hash(hash2_32(xi,   yi+1, seed), g01x, g01y);
    grad2_from_hash(hash2_32(xi+1, yi+1, seed), g11x, g11y);

    float d00 = g00x*tx       + g00y*ty;
    float d10 = g10x*(tx-1.f) + g10y*ty;
    float d01 = g01x*tx       + g01y*(ty-1.f);
    float d11 = g11x*(tx-1.f) + g11y*(ty-1.f);

    float ux = fade_f(tx), uy = fade_f(ty);
    float nx0 = lerp_f(d00, d10, ux);
    float nx1 = lerp_f(d01, d11, ux);
    float n   = lerp_f(nx0, nx1, uy);

    return 0.70710678f * n; // ~[-1,1]
}

// =================== fBm standard ===================
static inline float fbm_perlin(float x, float y, const HeightParams& P){
    float f = P.freq, a = P.amp, sum = 0.0f, norm = 0.0f;
    for(int o=0; o<P.octaves; ++o){
        sum  += a * perlin2D(x*f, y*f, P.seed + 131*o);
        norm += a;
        f *= P.lacunarity;
        a *= P.gain;
    }
    return (norm>1e-6f) ? (sum/norm) : 0.0f;
}

// =================== fBm "ridged" (crêtes) ===================
static inline float fbm_ridged(float x, float y, const HeightParams& P){
    float f = P.freq, a = P.amp, sum = 0.0f, norm = 0.0f;
    for(int o=0; o<P.octaves; ++o){
        float n = perlin2D(x*f, y*f, P.seed + 131*o);
        n = 1.0f - std::fabs(n);
        n *= n;
        sum  += a * n;
        norm += a;
        f *= P.lacunarity;
        a *= P.gain;
    }
    return (norm>1e-6f) ? (sum/norm) : 0.0f;
}

// =================== Domain warp léger ===================
static inline void warp_coords(float& x, float& y, const HeightParams& P){
    if (P.warpStrength <= 0.0f) return;
    float wx = perlin2D(x * P.warpFreq, y * P.warpFreq, P.seed+999);
    float wy = perlin2D((x+37.1f) * P.warpFreq, (y-11.9f) * P.warpFreq, P.seed-777);
    x += P.warpStrength * wx;
    y += P.warpStrength * wy;
}

static inline float hash01_u32(uint32_t h){
    h = (h ^ (h >> 13)) * 1274126177u;
    h ^= (h >> 16);
    return h * (1.0f/4294967295.0f);
}
static inline float hash2i01(int x, int y, int seed){
    uint32_t h = uint32_t(x)*374761393u + uint32_t(y)*668265263u + uint32_t(seed)*2246822519u;
    return hash01_u32(h);
}
static inline float valueNoise2D(float x, float y, int seed){
    int xi = int(std::floor(x)), yi = int(std::floor(y));
    float tx = x - xi, ty = y - yi;
    float a = hash2i01(xi,  yi,  seed);
    float b = hash2i01(xi+1,yi,  seed);
    float c = hash2i01(xi,  yi+1,seed);
    float d = hash2i01(xi+1,yi+1,seed);
    float ux = fade_f(tx), uy = fade_f(ty);
    float v = lerp_f( lerp_f(a,b,ux), lerp_f(c,d,ux), uy );
    return v*2.0f - 1.0f;
}
static inline float fbm_value(float x, float y, const HeightParams& P){
    float f=P.freq, a=P.amp, sum=0.f, norm=0.f;
    for(int o=0;o<P.octaves;++o){
        sum  += a * valueNoise2D(x*f, y*f, P.seed + 131*o);
        norm += a;
        f *= P.lacunarity;
        a *= P.gain;
    }
    return (norm>1e-6f) ? (sum/norm) : 0.0f;
}


// =================== Hauteur finale COSTAUD ===================
static inline float SampleTerrainHeight(float x, float y, const HeightParams& P){
    float X=x, Y=y;
    warp_coords(X,Y,P);

    float macro;
    if (P.model == NoiseType::Value) {
        HeightParams vP = P; vP.amp=1.0f; vP.freq=1.0f; vP.octaves=4; vP.lacunarity=2.0f; vP.gain=0.5f;
        macro = fbm_value(X*0.20f, Y*0.20f, vP);
    } else if (P.model == NoiseType::PerlinRidged) {
        HeightParams rP = P; rP.amp=1.0f; rP.freq=1.0f; rP.octaves=4; rP.lacunarity=2.0f; rP.gain=0.5f;
        float rid = fbm_ridged(X*0.20f, Y*0.20f, rP);
        macro = rid;
    } else { 
        HeightParams pP = P; pP.amp=1.0f; pP.freq=1.0f; pP.octaves=4; pP.lacunarity=2.0f; pP.gain=0.5f;
        macro = fbm_perlin(X*0.20f, Y*0.20f, pP);
    }

    float detail = (P.model==NoiseType::Value)
                 ? fbm_value (X*0.8f, Y*0.8f, P)
                 : fbm_perlin(X*0.8f, Y*0.8f, P);

    float basin = 0.0f;
    if (P.useBasins){
        basin = perlin2D(X * 0.28f, Y * 0.28f, P.seed + 777);
        basin = 0.5f*(basin+1.0f);
        basin = basin*basin*(3.0f - 2.0f*basin);
    }

    float h = 1.6f * macro + 0.22f * detail - (P.useBasins ? 2.2f * basin : 0.0f);

    return h * (P.amp / 10.0f) + P.heightBias;
}


