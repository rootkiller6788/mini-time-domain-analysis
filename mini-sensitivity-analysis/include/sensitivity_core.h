/**
 * sensitivity_core.h — Core Definitions for Sensitivity Analysis
 *
 * Domain: Time-Domain Analysis → Sensitivity Analysis (Control Theory)
 *
 * L1 Definitions:
 *   - Sensitivity function S(s) = 1 / (1 + L(s))
 *   - Complementary sensitivity T(s) = L(s) / (1 + L(s))
 *   - Loop transfer function L(s) = G(s) * H(s)
 *   - Return difference F(s) = 1 + L(s)
 *   - Absolute sensitivity S_abs
 *   - Relative sensitivity S_rel
 *   - Peak sensitivity Ms
 *   - Peak complementary sensitivity Mt
 *
 * L2 Core Concepts:
 *   - S + T = 1 (fundamental identity)
 *   - Waterbed effect: reducing sensitivity in one band increases it elsewhere
 *   - Disturbance rejection ∝ |S(jω)| small at low freq
 *   - Noise attenuation ∝ |T(jω)| small at high freq
 *   - Fundamental trade-off in feedback design
 *
 * References:
 *   - Bode, H.W. "Network Analysis and Feedback Amplifier Design" (1945)
 *   - Horowitz, I.M. "Synthesis of Feedback Systems" (1963)
 *   - Skogestad & Postlethwaite "Multivariable Feedback Control" (2005)
 */

#ifndef SENSITIVITY_CORE_H
#define SENSITIVITY_CORE_H

#include <stddef.h>
#include <math.h>

/* M_PI not defined in strict C11; define if missing */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Core Data Types — Complex Number Representation
 * ========================================================================== */

/**
 * Complex number in Cartesian form z = re + j*im.
 * Used for frequency-domain evaluation s = σ + jω.
 */
typedef struct {
    double re;   /**< Real part */
    double im;   /**< Imaginary part */
} Complex;

/**
 * Frequency response point: magnitude and phase at a given frequency.
 * |G(jω)| = mag, ∠G(jω) = phase (radians).
 */
typedef struct {
    double omega;     /**< Frequency in rad/s */
    double mag;       /**< Magnitude (linear, not dB) */
    double phase;     /**< Phase in radians */
} FreqPoint;

/* ==========================================================================
 * L1: Transfer Function Representation
 * ========================================================================== */

/** Maximum polynomial order for transfer function numerator/denominator */
#define MAX_TF_ORDER 20

/**
 * Rational transfer function: G(s) = b_n s^n + ... + b_1 s + b_0
 *                                    -------------------------------
 *                                    a_m s^m + ... + a_1 s + a_0
 * Stored as numerator coefficients b[0..n] and denominator a[0..m].
 * b[0] = constant term, b[n] = highest power coefficient.
 */
typedef struct {
    int n;                         /**< Degree of numerator polynomial */
    int m;                         /**< Degree of denominator polynomial */
    double b[MAX_TF_ORDER + 1];    /**< Numerator coefficients, b[0]=constant */
    double a[MAX_TF_ORDER + 1];    /**< Denominator coefficients, a[0]=constant */
} TransferFunction;

/* ==========================================================================
 * L1: Sensitivity Function Data
 * ========================================================================== */

/**
 * Sensitivity analysis result for a single frequency ω.
 * Contains S(jω), T(jω), and derived quantities.
 */
typedef struct {
    double omega;         /**< Frequency in rad/s */
    Complex S;            /**< Sensitivity function S(jω) = 1/(1+L(jω)) */
    Complex T;            /**< Complementary sensitivity T(jω) = L(jω)/(1+L(jω)) */
    double mag_S;         /**< Magnitude |S(jω)| */
    double mag_T;         /**< Magnitude |T(jω)| */
    double mag_S_dB;      /**< |S(jω)| in dB = 20·log10(|S|) */
    double mag_T_dB;      /**< |T(jω)| in dB = 20·log10(|T|) */
} SensitivityPoint;

/**
 * Full sensitivity analysis over a frequency sweep.
 */
typedef struct {
    int n_points;                      /**< Number of frequency points */
    SensitivityPoint *points;          /**< Array of points at each frequency */
    double Ms;                         /**< Peak sensitivity = max |S(jω)| */
    double Mt;                         /**< Peak complementary sensitivity = max |T(jω)| */
    double omega_Ms;                   /**< Frequency where Ms occurs */
    double omega_Mt;                   /**< Frequency where Mt occurs */
    double bandwidth;                  /**< -3dB bandwidth of T(s) (rad/s) */
} SensitivityAnalysis;

/* ==========================================================================
 * L1: Parameter Sensitivity Structures
 * ========================================================================== */

/**
 * Parameter sensitivity result: ∂y/∂p where y is output, p is parameter.
 * Measures how much output changes per unit parameter change.
 */
typedef struct {
    const char *param_name;      /**< Name of the parameter */
    double nominal_value;        /**< Nominal parameter value */
    double sensitivity_abs;      /**< Absolute sensitivity ∂y/∂p */
    double sensitivity_rel;      /**< Relative sensitivity (p/y)·(∂y/∂p) */
    double sensitivity_pct;      /**< Percentage sensitivity: relative × 100% */
} ParameterSensitivity;

/**
 * Collection of parameter sensitivities for a system.
 */
typedef struct {
    int n_params;                     /**< Number of parameters */
    ParameterSensitivity *params;     /**< Array of parameter sensitivities */
    double max_abs;                   /**< Maximum |absolute sensitivity| */
    int max_index;                    /**< Index of most sensitive parameter */
} ParamSensitivitySet;

/* ==========================================================================
 * L1: State-Space Representation (for time-domain sensitivity)
 * ========================================================================== */

/** Maximum state dimension */
#define MAX_STATE_DIM 20

/**
 * Linear time-invariant state-space model:
 *   ẋ = A x + B u
 *   y = C x + D u
 *
 * Matrices stored as flat row-major arrays for flexible index access.
 */
typedef struct {
    int n;           /**< Number of states */
    int m;           /**< Number of inputs */
    int p;           /**< Number of outputs */
    double A[MAX_STATE_DIM * MAX_STATE_DIM];  /**< State matrix (n×n), row-major */
    double B[MAX_STATE_DIM * MAX_STATE_DIM];  /**< Input matrix (n×m), row-major */
    double C[MAX_STATE_DIM * MAX_STATE_DIM];  /**< Output matrix (p×n), row-major */
    double D[MAX_STATE_DIM * MAX_STATE_DIM];  /**< Feedthrough matrix (p×m), row-major */
} StateSpace;

/* ==========================================================================
 * Complex Number Arithmetic (support for L3 mathematical structures)
 * ========================================================================== */

/**
 * Create complex number from real and imaginary parts.
 * Complexity: O(1)
 */
Complex complex_make(double re, double im);

/**
 * Complex addition: z1 + z2.
 * Complexity: O(1)
 */
Complex complex_add(Complex z1, Complex z2);

/**
 * Complex subtraction: z1 - z2.
 * Complexity: O(1)
 */
Complex complex_sub(Complex z1, Complex z2);

/**
 * Complex multiplication: z1 * z2.
 * Complexity: O(1)
 */
Complex complex_mul(Complex z1, Complex z2);

/**
 * Complex division: z1 / z2. Returns NaN on division by zero.
 * Complexity: O(1)
 */
Complex complex_div(Complex z1, Complex z2);

/**
 * Complex magnitude |z|.
 * Complexity: O(1)
 */
double complex_abs(Complex z);

/**
 * Complex argument (phase) in radians, range [-π, π].
 * Complexity: O(1)
 */
double complex_arg(Complex z);

/**
 * Complex conjugate: z̄ = re - j·im.
 * Complexity: O(1)
 */
Complex complex_conj(Complex z);

/**
 * Evaluate 1/(1+z) — the sensitivity function operator.
 * Complexity: O(1)
 */
Complex sensitivity_op(Complex z);

/**
 * Evaluate z/(1+z) — the complementary sensitivity operator.
 * Complexity: O(1)
 */
Complex comp_sensitivity_op(Complex z);

/* ==========================================================================
 * L2: Transfer Function Evaluation
 * ========================================================================== */

/**
 * Evaluate polynomial at complex s using Horner's method.
 * Exposed for use in checking pole/zero constraints.
 * p[0] is constant term, p[degree] is highest-power coefficient.
 * Complexity: O(degree)
 */
Complex poly_eval_horner(const double *coeff, int degree, Complex s);

/**
 * Evaluate transfer function G(s) at a complex frequency s = σ + jω.
 * Computes b(s)/a(s) using Horner's method for numerical stability.
 * Theorem: Horner's method (Horner, 1819) reduces evaluation from
 * O(n²) to O(n) and minimizes round-off error.
 *
 * Returns: G(s) as Complex. If denominator=0, returns (NaN,NaN).
 * Complexity: O(n + m) where n, m are polynomial degrees.
 */
Complex tf_evaluate(const TransferFunction *tf, Complex s);

/**
 * Evaluate loop transfer function L(s) = G(s) · H(s).
 * If H is NULL, treats H(s) = 1 (unity feedback).
 * Complexity: O(n_g + m_g + n_h + m_h)
 */
Complex loop_transfer(const TransferFunction *G, const TransferFunction *H, Complex s);

/**
 * Compute sensitivity S(s) = 1 / (1 + L(s)).
 * Handles the case L(s) → -1 (infinite sensitivity) by returning large value.
 * Complexity: O(1) after L(s) evaluation
 */
Complex compute_sensitivity(Complex L);

/**
 * Compute complementary sensitivity T(s) = L(s) / (1 + L(s)).
 * Complexity: O(1) after L(s) evaluation
 */
Complex compute_comp_sensitivity(Complex L);

/**
 * Verify S + T = 1 identity, returns |S+T-1|.
 * Should be ≤ 1e-12 for well-conditioned computations.
 * Complexity: O(1)
 */
double verify_st_identity(Complex S, Complex T);

/* ==========================================================================
 * Core Sensitivity Analysis API
 * ========================================================================== */

/**
 * Allocate and initialize a sensitivity analysis structure.
 * Caller must call sensitivity_analysis_free() when done.
 * Complexity: O(1) for initialization
 */
SensitivityAnalysis *sensitivity_analysis_alloc(int n_points);

/**
 * Free a sensitivity analysis structure and all associated memory.
 * Complexity: O(1)
 */
void sensitivity_analysis_free(SensitivityAnalysis *sa);

/**
 * Perform full sensitivity analysis of loop transfer L(s) over
 * logarithmically spaced frequencies from omega_min to omega_max.
 *
 * Computes S(jω), T(jω), Ms, Mt, bandwidth for each frequency.
 *
 * Theorem: S + T = 1 (Fundamental feedback identity, Black 1934)
 *
 * Complexity: O(N·(n+m)) where N = n_points
 */
void sensitivity_sweep(TransferFunction *G, TransferFunction *H,
                       double omega_min, double omega_max, int n_points,
                       SensitivityAnalysis *sa);

/**
 * Find the peak sensitivity Ms = max_ω |S(jω)| and the frequency at which
 * it occurs. Uses golden-section search for unimodal functions.
 * Complexity: O(log(N)) per search iteration
 */
double find_peak_sensitivity(TransferFunction *G, TransferFunction *H,
                             double omega_min, double omega_max,
                             double *omega_peak);

/**
 * Compute the -3dB bandwidth of the complementary sensitivity T(s).
 * Returns the frequency where |T(jω)| drops by 3dB from its DC value.
 * Uses bisection method with guaranteed convergence.
 * Complexity: O(log(ω_max/ω_min / ε)) iterations
 */
double compute_bandwidth(TransferFunction *G, TransferFunction *H,
                         double omega_max_search);

#ifdef __cplusplus
}
#endif

#endif /* SENSITIVITY_CORE_H */
