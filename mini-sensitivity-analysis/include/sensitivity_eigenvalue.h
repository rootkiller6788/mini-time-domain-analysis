/**
 * sensitivity_eigenvalue.h — Eigenvalue & Eigenvector Sensitivity Analysis
 *
 * Analyzes how eigenvalues λ_i and eigenvectors v_i of a matrix A
 * change with respect to parameter variations.
 *
 * L1 Definitions:
 *   - Eigenvalue sensitivity: ∂λ_k/∂p
 *   - Eigenvector sensitivity: ∂v_k/∂p
 *   - Eigenvalue condition number κ(λ_k) = 1/|y_k^* · x_k|
 *     where x_k is right eigenvector, y_k is left eigenvector
 *
 * L2 Core Concepts:
 *   - Eigenvalues determine system stability and time constants
 *   - Repeated eigenvalues have fundamentally different sensitivity behavior
 *   - Conditioning of eigenvalue problem affects sensitivity
 *
 * L4 Fundamental Laws (Theorems):
 *   - Wilkinson's formula: ∂λ/∂A_{ij} = (y_i · x_j) / (y^*·x)
 *     where x = right eigenvector, y = left eigenvector
 *   - Bauer-Fike theorem: |λ̂ - λ| ≤ κ_2(X)·‖ΔA‖_2
 *     for simple eigenvalue λ of A, where X is eigenvector matrix
 *
 *   - For a parameter p: ∂λ/∂p_i = (y^* · (∂A/∂p_i) · x) / (y^*·x)
 *     This follows from differentiating Ax = λx.
 *
 * L5 Computational Methods:
 *   - QR algorithm for eigenvalue decomposition
 *   - Power iteration for dominant eigenvalue
 *   - Inverse iteration for selected eigenvalues
 *   - Numerical differentiation of eigenvalues
 *
 * References:
 *   - Wilkinson, J.H. "The Algebraic Eigenvalue Problem" (1965)
 *   - Golub & Van Loan "Matrix Computations" (2013)
 *   - Stewart & Sun "Matrix Perturbation Theory" (1990)
 */

#ifndef SENSITIVITY_EIGENVALUE_H
#define SENSITIVITY_EIGENVALUE_H

#include "sensitivity_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Eigenvalue Sensitivity Structures
 * ========================================================================== */

/**
 * Eigenvalue with associated sensitivity information.
 */
typedef struct {
    Complex lambda;               /**< Eigenvalue λ = σ + jω */
    double sigma;                 /**< Real part (determines stability: σ<0 = stable) */
    double omega;                 /**< Imaginary part (determines oscillation freq) */
    double damping;               /**< Damping ratio ζ = -σ/|λ| (for complex eigenvalues) */
    double natural_freq;          /**< Natural frequency ω_n = |λ| */
    double *sensitivities;        /**< ∂λ/∂p_i for each parameter (n_params entries) */
    double *right_eigenvector;    /**< Right eigenvector x (n entries) */
    double *left_eigenvector;     /**< Left eigenvector y (n entries) */
    double condition_number;      /**< κ(λ) = ‖x‖·‖y‖/|y^*·x| (Wilkinson condition) */
} EigenvalueWithSensitivity;

/**
 * Full eigenvalue sensitivity analysis result.
 */
typedef struct {
    int n;                                    /**< Matrix dimension */
    int n_params;                             /**< Number of parameters */
    double **A_nominal;                       /**< Nominal matrix A (n×n) */
    double **dA_dp;                           /**< ∂A/∂p_i matrices (n_params × n×n), each flat */
    EigenvalueWithSensitivity *eigenvalues;    /**< n eigenvalues with sensitivities */
    double spectral_abscissa;                 /**< max Re(λ_i) — stability indicator */
    double spectral_radius;                   /**< max |λ_i| */
    double condition_number_A;                /**< 2-norm condition of eigenvector matrix */
} EigenSensitivityReport;

/* ==========================================================================
 * L5: Core Eigenvalue Sensitivity Computation
 * ========================================================================== */

/**
 * Compute eigenvalue sensitivity using Wilkinson's formula:
 *   ∂λ_k/∂p_i = (y_k^* · (∂A/∂p_i) · x_k) / (y_k^* · x_k)
 *
 * Requires the eigenvalue λ_k and its left (y_k) and right (x_k) eigenvectors,
 * plus the matrix derivative ∂A/∂p_i.
 *
 * Theorem: Wilkinson (1965) — eigenvalue derivative for simple eigenvalues.
 *
 * @param n matrix dimension
 * @param dA_dp partial derivative matrix ∂A/∂p (n×n, row-major)
 * @param lambda eigenvalue (Complex)
 * @param x right eigenvector (n entries)
 * @param y left eigenvector (n entries)
 * @return ∂λ/∂p as Complex number
 * Complexity: O(n²) for matrix-vector products
 */
Complex eigenvalue_param_sensitivity(int n, const double *dA_dp,
                                     Complex lambda,
                                     const double *x, const double *y);

/**
 * Compute eigenvalue sensitivity with respect to a specific matrix element A_{ij}:
 *   ∂λ/∂A_{ij} = (y_i · x_j) / (y^*·x)
 *
 * This is a special case where the derivative matrix has 1 at (i,j) and 0 elsewhere.
 * Complexity: O(1)
 */
double eigenvalue_matrix_sensitivity(int n, int i, int j,
                                     const double *x, const double *y);

/**
 * Compute all eigenvalue sensitivities for all parameters.
 * Fills the sensitivities array in EigenvalueWithSensitivity for each eigenvalue.
 *
 * @param report EigenSensitivityReport with eigenvalues and dA_dp populated
 * Complexity: O(n·n_params·n²)
 */
void compute_all_eigenvalue_sensitivities(EigenSensitivityReport *report);

/* ==========================================================================
 * L5: Eigenvalue Decomposition (QR Algorithm)
 * ========================================================================== */

/**
 * Compute all eigenvalues of a real n×n matrix using the QR algorithm
 * with Hessenberg reduction and double-shift strategy.
 *
 * Returns eigenvalues as an array of Complex numbers.
 * For real eigenvalues, im = 0.
 * For complex conjugate pairs, both eigenvalues are returned.
 *
 * Algorithm: Francis QR double-shift with Wilkinson shifts (Francis, 1961).
 *
 * @param A input matrix (n×n, row-major), modified in-place
 * @param n matrix dimension
 * @param eigenvalues output array of n Complex numbers
 * @return 0 on success, -1 if algorithm fails to converge
 * Complexity: O(n³)
 */
int qr_eigenvalues(double *A, int n, Complex *eigenvalues);

/**
 * Compute the eigenvectors of a real n×n matrix given its eigenvalues.
 * Uses inverse iteration with Rayleigh quotient refinement.
 *
 * @param A original matrix (n×n, row-major)
 * @param n dimension
 * @param eigenvalues previously computed eigenvalues
 * @param right_eigenvectors output (n×n, row-major: i-th row = i-th right eigenvector)
 * @param left_eigenvectors output (n×n, row-major: i-th row = i-th left eigenvector)
 * @return 0 on success, -1 on failure
 * Complexity: O(n³)
 */
int compute_eigenvectors(const double *A, int n, const Complex *eigenvalues,
                         double *right_eigenvectors, double *left_eigenvectors);

/**
 * Compute eigenvalue condition number (Wilkinson):
 *   κ(λ_k) = 1/|y_k^* · x_k|
 *
 * where x_k is the right eigenvector and y_k is the left eigenvector,
 * normalized so ‖x_k‖ = ‖y_k‖ = 1.
 *
 * Large κ indicates high sensitivity to perturbations.
 * Complexity: O(n)
 */
double eigenvalue_condition_number(int n, const double *x, const double *y);

/* ==========================================================================
 * L5: Power Iteration & Dominant Eigenvalues
 * ========================================================================== */

/**
 * Power iteration to find the dominant eigenvalue and eigenvector.
 * x_{k+1} = A·x_k / ‖A·x_k‖
 * λ ≈ (x_k^*·A·x_k) / (x_k^*·x_k) after convergence.
 *
 * Convergence rate: O(|λ_2/λ_1|^k) where |λ_1| > |λ_2| ≥ ...
 *
 * @param A matrix (n×n, row-major)
 * @param n dimension
 * @param max_iter maximum iterations
 * @param tol convergence tolerance
 * @param lambda output dominant eigenvalue
 * @param eigenvector output eigenvector (n entries, normalized)
 * @return number of iterations, or -1 if failed
 * Complexity: O(iterations·n²)
 */
int power_iteration(const double *A, int n, int max_iter, double tol,
                    double *lambda, double *eigenvector);

/**
 * Inverse iteration to find eigenvalue closest to a given shift μ.
 * (A - μI)^{-1}·x_{k+1} ∝ x_k
 *
 * @param A matrix (n×n, row-major)
 * @param n dimension
 * @param mu shift value
 * @param max_iter maximum iterations
 * @param tol convergence tolerance
 * @param lambda output eigenvalue nearest to mu
 * @param eigenvector output eigenvector (n entries, normalized)
 * @return iterations, or -1 if failed
 * Complexity: O(iterations·n³) for linear solve each step
 */
int inverse_iteration(const double *A, int n, double mu,
                      int max_iter, double tol,
                      double *lambda, double *eigenvector);

/* ==========================================================================
 * L5: Matrix Exponential for Time-Domain Sensitivity
 * ========================================================================== */

/**
 * Compute matrix exponential e^{A·t} using scaling-and-squaring with Padé
 * approximation. This is essential for time-domain sensitivity analysis
 * where the state transition matrix Φ(t) = e^{A·t} appears.
 *
 * Algorithm: Higham (2005) "The Scaling and Squaring Method..."
 *
 * @param A matrix (n×n, row-major)
 * @param n dimension
 * @param t time
 * @param expAt output matrix (n×n, row-major)
 * Complexity: O(n³)
 */
void matrix_exponential(const double *A, int n, double t, double *expAt);

/**
 * Compute the integral of matrix exponential:
 *   ∫_0^t e^{A·τ} dτ
 *
 * Useful for computing step response and its sensitivity:
 * y(t) = C·∫_0^t e^{A·τ}·B dτ · u0
 *
 * @param A matrix (n×n, row-major)
 * @param n dimension
 * @param t time
 * @param integral output matrix (n×n, row-major)
 * Complexity: O(n³)
 */
void matrix_exponential_integral(const double *A, int n, double t,
                                 double *integral);

/* ==========================================================================
 * L4: Stability Metrics from Eigenvalues
 * ========================================================================== */

/**
 * Compute spectral abscissa: α(A) = max_i Re(λ_i).
 * α < 0 ⇒ asymptotic stability for LTI systems.
 * Complexity: O(n)
 */
double spectral_abscissa(const Complex *eigenvalues, int n);

/**
 * Compute spectral radius: ρ(A) = max_i |λ_i|.
 * For discrete-time systems: ρ < 1 ⇒ stability.
 * Complexity: O(n)
 */
double spectral_radius(const Complex *eigenvalues, int n);

/**
 * Compute the damping ratio for each eigenvalue:
 * ζ_i = -Re(λ_i) / |λ_i| (for complex eigenvalues)
 * Returns array of n damping ratios.
 * Complexity: O(n)
 */
void compute_damping_ratios(const Complex *eigenvalues, int n,
                            double *damping);

/**
 * Estimate the settling time from eigenvalues:
 * t_s ≈ 4 / min_i |Re(λ_i)| for stable systems.
 * Based on the dominant (slowest) eigenvalue.
 * Complexity: O(n)
 */
double estimate_settling_time(const Complex *eigenvalues, int n);

/* ==========================================================================
 * Internal Helpers (exposed for testing)
 * ========================================================================== */

/**
 * Hessenberg reduction of a matrix. Transforms A to upper Hessenberg form
 * via Householder reflections. A = Q·H·Q^T.
 *
 * @param A input/output matrix (n×n, row-major)
 * @param n dimension
 * @param Q output orthogonal matrix (n×n, row-major)
 * Complexity: O(n³)
 */
void hessenberg_reduction(double *A, int n, double *Q);

/**
 * QR decomposition of a matrix using Householder reflections.
 * A = Q·R where Q is orthogonal, R is upper triangular.
 *
 * @param A input matrix (m×n, row-major), overwritten with R
 * @param m rows
 * @param n cols (m ≥ n)
 * @param Q output (m×m, row-major)
 * Complexity: O(m·n²)
 */
void qr_decomposition(double *A, int m, int n, double *Q);

#ifdef __cplusplus
}
#endif

#endif /* SENSITIVITY_EIGENVALUE_H */
