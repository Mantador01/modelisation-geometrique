#pragma once
#include "implicits.h"
#include "meshcolor.h"
#include "mathematics.h"

#include <vector>
#include <algorithm>
#include <cmath>


class BlobField : public AnalyticScalarField {
public:
    struct Primitive {
        enum Type { Point, Segment } type = Point;
        // Point
        Vector p;
        // Segment
        Vector a, b;
        // Influence
        double R = 0.25;
        double w = 1.0;

        Color  color{0.8,0.8,0.8};
    };

    void Clear() 
    {
        prims.clear();
    }

    void SetIso(double t) 
    {
        iso = t;
    }

    void AddPoint(const Vector& P, double R, double w, const Color& c = Color(0.8,0.8,0.8)) 
    {
        Primitive pr; pr.type = Primitive::Point;
        pr.p = P; pr.R = R; pr.w = w; pr.color = c;
        prims.push_back(pr);
    }

    void AddSegment(const Vector& A, const Vector& B, double R, double w,
                    const Color& c = Color(0.8,0.8,0.8)) {
        Primitive pr; pr.type = Primitive::Segment;
        pr.a = A; pr.b = B; pr.R = R; pr.w = w; pr.color = c;
        prims.push_back(pr);
    }

    double Value(const Vector& x) const override 
    {
        if (prims.empty())
            return Norm(x) - 1.0;

        double sum = 0.0;
        for (const auto& pr : prims) {
            double r = 0.0;
            if (pr.type == Primitive::Point) {
                r = Norm(x - pr.p);
            } else {
                Vector ab = pr.b - pr.a;
                double denom = Dot(ab, ab);
                double t = (denom > 0.0) ? Dot(x - pr.a, ab) / denom : 0.0;
#if __cplusplus >= 201703L
                t = std::clamp(t, 0.0, 1.0);
#else
                if (t < 0.0) t = 0.0; else if (t > 1.0) t = 1.0;
#endif
                Vector c = pr.a + t * ab;
                r = Norm(x - c);
            }
            if (r < pr.R) {
                double tt = 1.0 - (r*r) / (pr.R*pr.R);
                sum += pr.w * tt*tt*tt;
            }
        }
        return iso - sum;
    }

    const std::vector<Primitive>& Prims() const { return prims; }
    double Iso() const { return iso; }

private:
    std::vector<Primitive> prims;
    double iso = 0.5;
};
