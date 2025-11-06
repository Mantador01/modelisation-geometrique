#pragma once
#include "implicits.h"
#include <vector>

struct ErosionSphere {
    Vector c;
    double R;
    double smoothK;

    ErosionSphere(const Vector& C, double r, double k = 0.0)
      : c(C), R(r), smoothK(k) {}
};

class ErosionField : public AnalyticScalarField
{
public:
    const AnalyticScalarField* base = nullptr;
    std::vector<ErosionSphere> impacts;

    ErosionField() = default;
    explicit ErosionField(const AnalyticScalarField* baseField) : base(baseField) {}

    void SetBase(const AnalyticScalarField* baseField) { base = baseField; }
    void ClearImpacts() { impacts.clear(); }
    void AddImpact(const Vector& c, double R, double smoothK = 0.0) { impacts.push_back({c,R,smoothK}); }

    virtual double Value(const Vector& p) const override;

private:
    static inline double smin(double a, double b, double k);
    static inline double smax(double a, double b, double k);
};


// POUR TEST V2 : version avec mode érode/dépose

// #pragma once
// #include "implicits.h"
// #include <vector>

// enum class ImpactMode { Erode, Deposit };

// struct ErosionSphere {
//     Vector c;
//     double R;
//     double smoothK; 
//     ImpactMode mode; 

//     ErosionSphere(const Vector& C, double r, double k=0.0, ImpactMode m=ImpactMode::Erode)
//       : c(C), R(r), smoothK(k), mode(m) {}
// };

// class ErosionField : public AnalyticScalarField {
// public:
//     const AnalyticScalarField* base = nullptr;
//     std::vector<ErosionSphere> impacts;

//     ErosionField() = default;
//     explicit ErosionField(const AnalyticScalarField* baseField) : base(baseField) {}

//     void SetBase(const AnalyticScalarField* baseField) { base = baseField; }
//     void ClearImpacts() { impacts.clear(); }
//     void AddImpact(const Vector& c, double R, double k=0.0, ImpactMode m=ImpactMode::Erode)
//     { impacts.push_back({c,R,k,m}); }

//     double Value(const Vector& p) const override;

// private:
//     static inline double smin(double a, double b, double k) {
//         if (k<=0) return std::min(a,b);
//         double t = 0.5 + 0.5*(b - a)/k;
//         double h = (t < 0.0) ? 0.0 : (t > 1.0) ? 1.0 : t;
//         return (b*(1.0 - h) + a*h) - k*h*(1.0 - h);
//     }
//     static inline double smax(double a, double b, double k) {
//         return -smin(-a, -b, k);
//     }
// };
