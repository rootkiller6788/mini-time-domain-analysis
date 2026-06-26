#ifndef SENSITIVITY_ANALYSIS_H
#define SENSITIVITY_ANALYSIS_H

#include "steady_state_error.h"

/* ============================================================================
 * mini-steady-state-error: Sensitivity Analysis Module
 *
 * L2: Sensitivity function definition: S(s) = 1/(1 + G(s)H(s))
 * L3: Relationship between sensitivity and steady-state error
 * L5: Parameter sensitivity, robustness margins
 * L6: Bode integral (waterbed effect)
 * L8: Structured singular value (mu) for robust steady-state error
 *
 * Reference: Nise §7.7, Ogata §5.8, Skogestad & Postlethwaite (2005)
 * ============================================================================ */

/* L2: Sensitivity function values at key frequencies */
typedef struct {
    double S_dc;           /* |S(0)| = sensitivity at DC */
    double S_peak;         /* max |S(jω)| = M_s (robustness indicator) */
    double S_bandwidth;    /* frequency where |S| crosses -3 dB */
    double T_dc;           /* |T(0)| = complementary sensitivity at DC */
    double T_peak;         /* max |T(jω)| = M_t */
} SensitivityMetrics;

/* L5: Parameter sensitivity of steady-state error */
typedef struct {
    double d_ess_d_K;      /* ∂e_ss/∂K : gain sensitivity */
    double d_ess_d_pole;   /* ∂e_ss/∂p : dominant pole sensitivity */
    double d_ess_d_zero;   /* ∂e_ss/∂z : zero location sensitivity */
    double worst_case_e_ss;/* e_ss under worst-case parameter variation */
    double param_variance; /* normalized variance-based sensitivity */
} ParameterSensitivity;

/* L6: Bode integral result */
typedef struct {
    double waterbed_integral;   /* integral of ln|S(jω)| dω */
    double positive_area;       /* area where |S| > 1 (0 dB) */
    double negative_area;       /* area where |S| < 1 (0 dB) */
    bool   waterbed_violated;   /* true if non-minimum-phase or unstable */
} BodeIntegral;

/* L8: Structured singular value analysis */
typedef struct {
    double mu_upper_bound;      /* upper bound of mu */
    double mu_lower_bound;      /* lower bound of mu */
    double robust_margin;       /* 1/mu: robust stability margin */
    double nominal_e_ss;        /* nominal steady-state error */
    double worst_case_e_ss_mu;  /* worst-case SSE under structured uncertainty */
    int    perturbation_count;  /* number of uncertainty blocks */
} StructuredSSV;

/* ---- L2: Sensitivity function computation ---- */

/**
 * L2: Compute sensitivity function S(s) = 1/(1 + L(s)) where L(s) = G(s)H(s).
 * S(s) represents the transfer function from reference to error.
 *
 * @param open_loop  open-loop transfer function L(s) = G(s)H(s)
 * @return           S(s) = 1/(1+L(s)) as TransferFunction (caller frees)
 *
 * Complexity: O((n+m)^2)
 */
TransferFunction sa_compute_sensitivity(const TransferFunction *open_loop);

/**
 * L2: Compute complementary sensitivity T(s) = L(s)/(1+L(s)).
 * T(s) = 1 - S(s) (for unity feedback).
 *
 * @param open_loop  open-loop transfer function L(s)
 * @return           T(s) as TransferFunction (caller frees)
 *
 * Complexity: O((n+m)^2)
 */
TransferFunction sa_compute_complementary_sensitivity(const TransferFunction *open_loop);

/**
 * L2: Compute DC sensitivity |S(0)|.
 * For Type 0: |S(0)| = 1/(1+Kp)
 * For Type ≥1: |S(0)| = 0 (integral action forces zero DC error)
 *
 * Complexity: O(1)
 */
double sa_dc_sensitivity(const ErrorConstants *ec);

/**
 * L3: Compute sensitivity at a specific frequency ω.
 * |S(jω)| = 1 / |1 + L(jω)|
 *
 * @param open_loop  open-loop transfer function
 * @param omega      frequency (rad/s)
 * @return           |S(jω)|
 *
 * Complexity: O(n+m) for polynomial evaluation at jω
 */
double sa_sensitivity_at_freq(const TransferFunction *open_loop, double omega);

/* ---- L5: Parameter sensitivity ---- */

/**
 * L5: Compute sensitivity of steady-state error to gain variations.
 * d(e_ss)/dK for a Type 0 system with step input:
 *   e_ss = 1/(1+K) → de_ss/dK = -1/(1+K)^2 = -e_ss^2
 *
 * @param system_type  system type
 * @param input_type   test input type
 * @param K            current gain value
 * @param e_ss_current current steady-state error
 * @return             d(e_ss)/dK
 *
 * Complexity: O(1)
 */
double sa_gain_sensitivity(int system_type, TestInputType input_type,
                            double K, double e_ss_current);

/**
 * L5: Compute parameter sensitivity vector for all dominant parameters.
 * Uses finite difference approximation to compute ∂e_ss/∂p_i for each parameter.
 *
 * @param G_nominal     nominal plant transfer function
 * @param param_delta   fractional parameter perturbation (e.g., 0.01 = 1%)
 * @return              filled ParameterSensitivity struct
 *
 * Complexity: O(k * n^2) where k = number of parameters
 */
ParameterSensitivity sa_parameter_sensitivity(const TransferFunction *G_nominal,
                                                double param_delta);

/**
 * L5: Compute condition number for error constant computation.
 * High condition number → e_ss is highly sensitive to small parameter changes.
 * κ = |d(e_ss)/d(param)| * |param / e_ss|
 *
 * Complexity: O(1)
 */
double sa_error_condition_number(double e_ss, double d_ess_d_param,
                                  double param_value);

/* ---- L6: Bode integral (waterbed effect) ---- */

/**
 * L6: Bode sensitivity integral theorem.
 * For stable open-loop L(s) with relative degree ≥ 2:
 *   ∫_0^∞ ln|S(jω)| dω = 0
 *
 * This means sensitivity reduction in one frequency range forces
 * sensitivity increase in another (the "waterbed effect").
 *
 * For unstable plants with RHP poles p_i:
 *   ∫_0^∞ ln|S(jω)| dω = π * Σ Re(p_i)
 *
 * This function approximates the integral using numerical quadrature
 * over a finite frequency range.
 *
 * @param open_loop    open-loop transfer function L(s)
 * @param rhp_poles    array of RHP pole locations (real parts > 0)
 * @param num_rhp      number of RHP poles
 * @param omega_min    lower integration bound (rad/s)
 * @param omega_max    upper integration bound (rad/s)
 * @param num_points   number of frequency points for quadrature
 * @return             filled BodeIntegral struct
 *
 * Reference: Bode (1945), Freudenberg & Looze (1985)
 * Complexity: O(N * (n+m)) for N frequency evaluations
 */
BodeIntegral sa_bode_integral(const TransferFunction *open_loop,
                               const double *rhp_poles, int num_rhp,
                               double omega_min, double omega_max,
                               int num_points);

/**
 * L6: Compute waterbed constraint for a given frequency range.
 * If |S| < 1 from 0 to ω_a (disturbance rejection band), then
 * there must exist ω > ω_a where |S| > 1 to satisfy Bode integral.
 *
 * @param reduction_db      sensitivity reduction in dB (negative = suppression)
 * @param bandwidth         bandwidth of reduction ω_a (rad/s)
 * @param relative_degree   relative degree of L(s)
 * @return                  minimum peak |S| required above ω_a
 *
 * Complexity: O(1) — analytical bound
 */
double sa_waterbed_peak_bound(double reduction_db, double bandwidth,
                               int relative_degree);

/* ---- L8: Structured singular value (mu) analysis ---- */

/**
 * L8: Compute structured singular value for robust steady-state error.
 *
 * μ(M) = 1 / min{σ̄(Δ) : det(I - MΔ) = 0, Δ ∈ Δ_structure}
 *
 * This quantifies the smallest structured uncertainty that can cause
 * the steady-state error to become unbounded or exceed specifications.
 *
 * Uses the Osborne iteration for μ upper bound computation.
 *
 * @param M                  M matrix (sensitivity of error to perturbations)
 * @param dim                dimension of M matrix (square)
 * @param block_structure    array of block sizes for structured Δ
 * @param num_blocks         number of uncertainty blocks
 * @param max_iter           maximum Osborne iterations
 * @param tol                convergence tolerance
 * @return                   filled StructuredSSV struct
 *
 * Reference: Doyle (1982), Packard & Doyle (1993)
 * Complexity: O(max_iter * dim^3)
 */
StructuredSSV sa_structured_sv(const double *M, int dim,
                                const int *block_structure, int num_blocks,
                                int max_iter, double tol);

/**
 * L8: Compute robust steady-state error margin.
 * Given nominal plant with multiplicative uncertainty Δ(s), bounded as |Δ(jω)| ≤ |W(jω)|,
 * compute the maximum possible steady-state error under worst-case uncertainty.
 *
 * @param G_nominal        nominal plant
 * @param W_uncertainty    uncertainty weight |W(jω)| bounds |Δ(jω)|
 * @param omega_range      frequency range of interest (2-element: [min, max])
 * @return                 worst-case steady-state error upper bound
 *
 * Complexity: O(N * n^2) for N frequency grid points
 */
double sa_robust_sse_bound(const TransferFunction *G_nominal,
                            const TransferFunction *W_uncertainty,
                            const double omega_range[2]);

/* ---- Utility ---- */

void sa_metrics_print(const SensitivityMetrics *m);
void sa_paramsens_print(const ParameterSensitivity *ps);
void sa_bode_print(const BodeIntegral *bi);

#endif /* SENSITIVITY_ANALYSIS_H */
