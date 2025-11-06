#pragma once
#include "implicits.h"
#include <memory>

struct Impact { Vector c; double R; double k; };

class DiffSphereField : public AnalyticScalarField {
public:
  const AnalyticScalarField* prev = nullptr;
  Impact imp;
  DiffSphereField(const AnalyticScalarField* P, Impact I) : prev(P), imp(I) {}

  static inline double smin(double a, double b, double k){
    double t = 0.5 + 0.5*(b - a)/k;
    double h = (t < 0.0) ? 0.0 : (t > 1.0) ? 1.0 : t;
    return (b*(1.0 - h) + a*h) - k*h*(1.0 - h);
  }
  static inline double smax(double a, double b, double k){ return -smin(-a, -b, k); }

  double Value(const Vector& p) const override {
    double v  = prev ? prev->Value(p) : (Norm(p) - 1.0);
    double sd = Norm(p - imp.c) - imp.R;
    double cut = -sd;
    return (imp.k>0.0)? smax(v, cut, imp.k) : std::max(v, cut);
  }
};


class ErosionIncremental {
public:
  std::unique_ptr<AnalyticScalarField> ownedBase; 
  const AnalyticScalarField* current = nullptr;

  explicit ErosionIncremental(const AnalyticScalarField* base) : current(base) {}
  void Own(AnalyticScalarField* base){ ownedBase.reset(base); current = ownedBase.get(); }

  void Apply(const Vector& c, double R, double k=0.0) {
    current = new DiffSphereField(current, {c,R,k});
  }
};
