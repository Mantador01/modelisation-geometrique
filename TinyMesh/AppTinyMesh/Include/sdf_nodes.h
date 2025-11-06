#pragma once
#include "blobfield.h"
#include "implicits.h"
#include <memory>
#include <algorithm> 

static inline double clamp01(double x){ return x<0.0?0.0:(x>1.0?1.0:x); }
static inline double lerp(double a, double b, double t){ return a + t*(b - a); }

static inline Vector absComp(const Vector& v){
    return Vector(std::abs(v[0]), std::abs(v[1]), std::abs(v[2]));
}
static inline Vector maxComp(const Vector& a, const Vector& b){
    return Vector(std::max(a[0], b[0]), std::max(a[1], b[1]), std::max(a[2], b[2]));
}

static inline double smin_poly(double a, double b, double k){
    const double h = clamp01(0.5 + 0.5*(b - a)/k);
    return lerp(b, a, h) - k*h*(1.0 - h);
}

static inline double lerp_(double a, double b, double t){ return a + t*(b-a); }
static inline double smin_poly_ncxx20(double a, double b, double k) {
    double h = clamp01(0.5 + 0.5*(b - a)/k);
    return lerp_(b, a, h) - k*h*(1.0 - h);
}


class SDFNode : public AnalyticScalarField {
public:
    virtual ~SDFNode() {}
    virtual double Value(const Vector& p) const = 0;
};

// ================= PRIMITIVES =================

class SphereNode : public SDFNode {
public:
    Vector center;
    double radius;

    SphereNode(const Vector& c, double r) : center(c), radius(r) {}
    double Value(const Vector& p) const override {
        return Norm(p - center) - radius;
    }
};

class BoxNode : public SDFNode {
public:
    Vector center, size; // half-extents
    BoxNode(const Vector& c, const Vector& s) : center(c), size(s) {}
    double Value(const Vector& p) const override {
        const Vector q = absComp(p - center) - size;
        const Vector qpos = maxComp(q, Vector(0.0));
        // outside distance + inside max component
        return Norm(qpos) + std::min(std::max(q[0], std::max(q[1], q[2])), 0.0);
    }
};

class CapsuleNode : public SDFNode {
public:
    Vector a, b;
    double r;
    CapsuleNode(const Vector& A, const Vector& B, double R) : a(A), b(B), r(R) {}
    double Value(const Vector& p) const override {
        const Vector pa = p - a, ba = b - a;
        const double denom = Dot(ba, ba);
        double h = (denom > 0.0) ? Dot(pa, ba) / denom : 0.0;
        h = clamp01(h);
        return Norm(pa - ba * h) - r;
    }
};

// Tore centré, axe Z, rayon majeur R (au centre du tube) et rayon mineur r (rayon du tube)
class TorusNode : public SDFNode {
public:
    Vector center;
    double R; // grand rayon
    double r; // petit rayon
    TorusNode(const Vector& c, double Rmajor, double Rminor) : center(c), R(Rmajor), r(Rminor) {}
    double Value(const Vector& p) const override {
        Vector q = p - center;
        // distance dans le plan XY au cercle de rayon R
        double qxy = std::sqrt(q[0]*q[0] + q[1]*q[1]) - R;
        return std::sqrt(qxy*qxy + q[2]*q[2]) - r;
    }
};


// ================= OPERATIONS =================

class UnionNode : public SDFNode {
public:
    std::shared_ptr<SDFNode> A, B;
    UnionNode(std::shared_ptr<SDFNode> a, std::shared_ptr<SDFNode> b) : A(a), B(b) {}
    double Value(const Vector& p) const override {
        return std::min(A->Value(p), B->Value(p));
    }
};

class IntersectionNode : public SDFNode {
public:
    std::shared_ptr<SDFNode> A, B;
    IntersectionNode(std::shared_ptr<SDFNode> a, std::shared_ptr<SDFNode> b) : A(a), B(b) {}
    double Value(const Vector& p) const override {
        return std::max(A->Value(p), B->Value(p));
    }
};

class DifferenceNode : public SDFNode {
public:
    std::shared_ptr<SDFNode> A, B;
    DifferenceNode(std::shared_ptr<SDFNode> a, std::shared_ptr<SDFNode> b) : A(a), B(b) {}
    double Value(const Vector& p) const override {
        return std::max(A->Value(p), -B->Value(p));
    }
};

class BlendNode : public SDFNode {
public:
    std::shared_ptr<SDFNode> A, B;
    double k;
    BlendNode(std::shared_ptr<SDFNode> a, std::shared_ptr<SDFNode> b, double K) : A(a), B(b), k(K) {}
    double Value(const Vector& p) const override {
        const double da = A->Value(p);
        const double db = B->Value(p);
        const double h = clamp01(0.5 + 0.5 * (db - da) / k);
        return lerp(db, da, h) - k * h * (1.0 - h);
    }
};

class BlobAsSDFNode : public SDFNode {
public:
    const BlobField* blob = nullptr; 
    explicit BlobAsSDFNode(const BlobField* b) : blob(b) {}
    double Value(const Vector& p) const override {
        return blob ? blob->Value(p) : 1e9;
    }
};

class BlobApproxSDFNode : public SDFNode {
public:
    const BlobField* blob = nullptr;
    double k; // largeur de lissage (ex: 0.15 * moyenne des R)

    BlobApproxSDFNode(const BlobField* b, double smoothK) : blob(b), k(smoothK) {}

    double Value(const Vector& x) const override {
        if (!blob) return 1e9;
        const auto& P = blob->Prims();
        if (P.empty()) return Norm(x) - 1.0;

        double d = 1e9;

        for (const auto& pr : P) {
            double di = 1e9;
            if (pr.type == BlobField::Primitive::Point) {
                di = Norm(x - pr.p) - pr.R; // sphère
            } else { // Segment -> capsule
                Vector ab = pr.b - pr.a;
                double denom = Dot(ab, ab);
                double t = (denom > 0.0) ? Dot(x - pr.a, ab) / denom : 0.0;
                t = clamp01(t);
                Vector c = pr.a + t * ab;
                di = Norm(x - c) - pr.R;
            }

#if __cplusplus >= 202002L
            d = smin_poly(d, di, k);
#else
            d = smin_poly_ncxx20(d, di, k);
#endif
        }
        return d;
    }
};
