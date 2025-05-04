This repository contains solutions to homework assignments for the Numerical Analysis and Computing course (RU). Each homework is organized with source code and brief explanations below.

HW1: Root‑Finding Methods
Implements three classic algorithms to find a root of a given function, stopping when successive approximations differ by less than 10⁻⁶.

Bisection Method halves an interval that brackets the root.
Newton’s Method uses the function and its derivative to converge quickly.
Secant Method approximates the derivative using two previous estimates, requiring no explicit derivative.

HW3: Newton Interpolation of eˣ
Builds Newton’s interpolating polynomial for eˣ on equally spaced nodes over [–1, 1]. Samples 501 points in the interval to measure the maximum difference between the true exponential and the polynomial approximation, reporting how the error changes as the polynomial degree increases.

HW4: LU Factorization Solver
Performs LU decomposition of a square matrix without pivoting, then uses forward and backward substitution to solve Ax = b. Reads a matrix and right‑hand vector and outputs the solution.

All code is in Python 3. 
