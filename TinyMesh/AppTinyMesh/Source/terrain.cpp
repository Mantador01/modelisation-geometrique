#include "terrain.h"
#include <vector>
#include <algorithm>

void BuildTerrainMesh(const HeightParams& P, Mesh& out)
{
    out = Mesh();

    const int N = P.N;
    const float s = P.scale;
    auto idx = [N](int x, int y) { return y * N + x; };

    std::vector<float> H(N * N, 0.0f);

    float x0 = 0.0f, z0 = 0.0f;
    if (P.center) {
        x0 = -(N - 1) * 0.5f * s;
        z0 = -(N - 1) * 0.5f * s;
    }

    for (int y = 0; y < N; ++y){
        for (int x = 0; x < N; ++x){
            float X = x0 + x * s;
            float Z = z0 + y * s;
            H[idx(x, y)] = SampleTerrainHeight(X, Z, P);
        }
    }

    std::vector<Vector> verts; verts.reserve(N*N);
    for (int y = 0; y < N; ++y){
        for (int x = 0; x < N; ++x){
            float X = x0 + x * s;
            float Z = z0 + y * s;
            float Y = H[idx(x, y)];
            verts.emplace_back(X, Y, Z);
        }
    }

    std::vector<int> indices; indices.reserve((N-1)*(N-1)*6);
    for (int y = 0; y < N - 1; ++y){
        for (int x = 0; x < N - 1; ++x){
            int i0 = idx(x, y), i1 = idx(x + 1, y);
            int i2 = idx(x, y + 1), i3 = idx(x + 1, y + 1);
            indices.push_back(i0); indices.push_back(i2); indices.push_back(i1);
            indices.push_back(i1); indices.push_back(i2); indices.push_back(i3);
        }
    }

    Mesh m(verts, indices);
    m.SmoothNormals();
    out = m;
}
