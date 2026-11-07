#pragma once
#include <cstddef>
#if defined(_WIN32)
# if defined(PC_BUILD)
#  define PC_API __declspec(dllexport)
# else
#  define PC_API __declspec(dllimport)
# endif
#else
# define PC_API __attribute__((visibility("default")))
#endif

// ABI 1: caller-owned nonoverlapping arrays, n entries except off[n-1].
// No pointers retained, no exceptions cross ABI, outputs untouched on failure.
// Status: 0 success, 1 invalid input, 2 unrepresentable arithmetic,
// 3 non-positive LDL pivot, 4 allocation/internal failure.
extern "C" {
PC_API int pc_abi_version() noexcept;
// capacity=rho*c, latent=rho*L. Residual [J], diagonal/off [J/K].
// side=0: right derivative at endpoints; side=1: left derivative.
// Negative trial T is allowed for line search; old T must be nonnegative.
PC_API int pc_assemble(std::size_t n, double capacity, double latent,
    double solidus, double liquidus, double conductivity, double dx,
    double area, double dt, int has_left, double left, int has_right,
    double right, int side, const double* old_t, const double* trial_t,
    const double* source, double* residual, double* diagonal, double* off) noexcept;
// Symmetric positive definite tridiagonal LDL^T solve. No pivoting/fallback.
PC_API int pc_solve(std::size_t n, const double* diagonal,
    const double* off, const double* rhs, double* solution) noexcept;
}
