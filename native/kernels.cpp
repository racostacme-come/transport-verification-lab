#include "kernels.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
bool positive(double x) { return std::isfinite(x) && x > 0.0; }
double fraction(double t, double ts, double width) {
    return std::clamp((t-ts)/width, 0.0, 1.0);
}
}
extern "C" int pc_abi_version() noexcept { return 1; }

extern "C" int pc_assemble(std::size_t n, double capacity, double latent,
    double ts, double tl, double conductivity, double dx, double area, double dt,
    int has_left, double left, int has_right, double right, int side,
    const double* old_t, const double* t, const double* source,
    double* residual, double* diagonal, double* off) noexcept {
    if (n == 0 || !old_t || !t || !source || !residual || !diagonal || (n > 1 && !off)
        || !positive(capacity) || !std::isfinite(latent) || latent < 0.0
        || !std::isfinite(ts) || ts < 0.0 || !std::isfinite(tl) || tl <= ts
        || !positive(conductivity) || !positive(dx) || !positive(area) || !positive(dt)
        || (side != 0 && side != 1) || (has_left != 0 && has_left != 1)
        || (has_right != 0 && has_right != 1)
        || (has_left && (!std::isfinite(left) || left < 0.0))
        || (has_right && (!std::isfinite(right) || right < 0.0))) return 1;
    const double width = tl-ts, volume = area*dx;
    const double face = dt*conductivity*area/dx, mushy = capacity+latent/width;
    if (!positive(volume) || !positive(face) || !positive(mushy)) return 2;
    for (std::size_t i=0; i<n; ++i)
        if (!std::isfinite(old_t[i]) || old_t[i]<0.0 || !std::isfinite(t[i])
            || !std::isfinite(source[i])) return 1;
    try {
        std::vector<double> r(n), d(n), e(n-1, -face);
        for (std::size_t i=0; i<n; ++i) {
            const double change = capacity*(t[i]-old_t[i])
                + latent*(fraction(t[i],ts,width)-fraction(old_t[i],ts,width));
            r[i] = volume*(change-dt*source[i]);
            const bool inside = side==0 ? (t[i]>=ts && t[i]<tl) : (t[i]>ts && t[i]<=tl);
            d[i] = volume*(inside ? mushy : capacity);
        }
        for (std::size_t i=0; i+1<n; ++i) {
            const double exchange = face*(t[i+1]-t[i]);
            r[i] -= exchange; r[i+1] += exchange;
            d[i] += face; d[i+1] += face;
        }
        if (has_left) { r[0] += 2*face*(t[0]-left); d[0] += 2*face; }
        if (has_right) { r[n-1] += 2*face*(t[n-1]-right); d[n-1] += 2*face; }
        for (std::size_t i=0; i<n; ++i)
            if (!std::isfinite(r[i]) || !positive(d[i])) return 2;
        std::copy(r.begin(),r.end(),residual); std::copy(d.begin(),d.end(),diagonal);
        if (n>1) std::copy(e.begin(),e.end(),off);
        return 0;
    } catch (...) { return 4; }
}

extern "C" int pc_solve(std::size_t n, const double* diagonal,
    const double* off, const double* rhs, double* solution) noexcept {
    if (n==0 || !diagonal || !rhs || !solution || (n>1 && !off)) return 1;
    for (std::size_t i=0; i<n; ++i)
        if (!std::isfinite(diagonal[i]) || !std::isfinite(rhs[i])
            || (i+1<n && !std::isfinite(off[i]))) return 1;
    try {
        std::vector<double> pivots(n), factors(n-1), x(rhs,rhs+n);
        pivots[0]=diagonal[0];
        if (pivots[0]<=0) return 3;
        for (std::size_t i=1; i<n; ++i) {
            factors[i-1]=off[i-1]/pivots[i-1];
            pivots[i]=diagonal[i]-factors[i-1]*off[i-1];
            x[i]-=factors[i-1]*x[i-1];
            if (!std::isfinite(factors[i-1]) || !std::isfinite(pivots[i]) || !std::isfinite(x[i])) return 2;
            if (pivots[i]<=0) return 3;
        }
        x[n-1]/=pivots[n-1];
        for (std::size_t i=n-1; i>0; --i) x[i-1]=x[i-1]/pivots[i-1]-factors[i-1]*x[i];
        for (double value:x) if (!std::isfinite(value)) return 2;
        std::copy(x.begin(),x.end(),solution);
        return 0;
    } catch (...) { return 4; }
}
