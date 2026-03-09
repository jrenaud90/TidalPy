#pragma once

#include <complex>

/// LAPACK zgesv: Solve a complex double-precision linear system A * X = B.
///
/// This is the standard Fortran calling convention for LAPACK's zgesv.
/// The symbol `zgesv_` (with trailing underscore) is provided by the system LAPACK
/// or by scipy's bundled OpenBLAS/MKL at link time.
///
/// Important notes:
///   - Matrix A must be in Fortran (column-major) order.
///   - A is overwritten with its LU factorization on output.
///   - B is overwritten with the solution X on output.
///   - ipiv is a pivot index array of size >= n.
///   - info == 0 on success, > 0 if singular, < 0 if bad argument.
extern "C" {
    void zgesv_(
        int* n,                        // (Input)  Order of matrix A (N x N)
        int* nrhs,                     // (Input)  Number of right-hand-side columns in B
        std::complex<double>* a,       // (In/Out) Coefficient matrix A, overwritten with LU
        int* lda,                      // (Input)  Leading dimension of A (>= n)
        int* ipiv,                     // (Output) Pivot indices of size >= n
        std::complex<double>* b,       // (In/Out) RHS matrix B, overwritten with solution X
        int* ldb,                      // (Input)  Leading dimension of B (>= n)
        int* info                      // (Output) 0=success, >0=singular, <0=bad arg
    );
}
