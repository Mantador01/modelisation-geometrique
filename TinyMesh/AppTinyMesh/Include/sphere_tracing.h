#pragma once
#include "implicits.h"
#include "ray.h"

bool IntersectSphereTracing(const AnalyticScalarField& sdf,
                            const Ray& ray,
                            double& tHit,
                            Vector* hitP = nullptr,
                            Vector* hitN = nullptr);
