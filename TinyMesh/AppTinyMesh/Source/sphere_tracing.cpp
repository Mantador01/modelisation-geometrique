#include "sphere_tracing.h"
#include <algorithm>
#include <cmath>

bool IntersectSphereTracing(const AnalyticScalarField& sdf,
                            const Ray& ray,
                            double& tHit,
                            Vector* hitP,
                            Vector* hitN)
{
    const double tMin     = 0.0;
    const double tMax     = 200.0;
    const int    maxSteps = 256;
    const double epsHit   = 1e-4;
    const double epsStep  = 1e-6;

    double t = tMin;
    for(int i=0; i<maxSteps && t < tMax; ++i)
    {
        Vector p = ray(t);
        double d = sdf.Value(p); 
        if(d < epsHit)
        {
            tHit = t;
            if(hitP) *hitP = p;
            if(hitN) *hitN = sdf.Normal(p); 
            return true;
        }
        t += std::max(d, epsStep);
    }
    return false;
}
