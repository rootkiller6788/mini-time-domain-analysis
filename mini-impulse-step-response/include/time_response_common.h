/**
 * @file time_response_common.h
 * @brief Common type definitions, system models, and response structures
 *        for time-domain impulse/step response analysis.
 *
 * L1 --- Definitions: System model structs (FirstOrderModel, SecondOrderModel,
 *       TransferFunction, StateSpaceModel), response containers
 * L2 --- Core Concepts: LTI system representation, transfer function <-> state-space,
 *       superposition, causality, BIBO stability in time domain
 * L3 --- Engineering Quantities: Time vectors, sample rates, tolerance thresholds
 * L4 --- Fundamental Laws: Convolution integral, Initial/Final Value Theorems
 *
 * Course alignment:
 *   MIT 6.302 Feedback Systems -- time-domain response characterization
 *   Stanford ENGR105 Feedback Control -- transient specifications
 *   Berkeley ME132 Dynamic Systems -- first/second-order time response
 *   Caltech CDS 110 Intro to Control -- impulse/step fundamentals
 *   ETH 151-0591 Control I -- Zeitantwortanalyse
 *   Cambridge 3F2 Systems & Control -- transient response
 *   Tsinghua Automatic Control Theory -- time-domain analysis method
 *
 * Textbook: Ogata, "Modern Control Engineering" (2010)
 *           Dorf & Bishop, "Modern Control Systems" (2017)
 *           Franklin, Powell, Emami-Naeini, "Feedback Control" (2019)
 */

#ifndef TIME_RESPONSE_COMMON_H
#define TIME_RESPONSE_COMMON_H

#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L1: Core System Model Types
 * ========================================================================== */

/**
 * First-order system model: G(s) = K / (tau*s + 1)
 *
 * K   -- steady-state gain (dimensionless or [output]/[input])
 * tau -- time constant [seconds]
 *
 * Step response:  y(t) = K * (1 - exp(-t/tau))
 * Impulse response: y(t) = (K/tau) * exp(-t/tau)
 *
 * The time constant tau is the time required to reach 63.2% of steady-state
 * for a step input. After 4*tau, the response reaches ~98.2% of final value.
 */
typedef struct {
    double K;      /**< DC gain (steady-state gain) */
    double tau;    /**< time constant [s] */
} FirstOrderModel;

/**
 * Second-order system model: G(s) = K * omega_n^2 / (s^2 + 2*zeta*omega_n*s + omega_n^2)
 *
 * K       -- DC gain
 * zeta    -- damping ratio (dimensionless)
 * omega_n -- natural frequency [rad/s]
 *
 * Pole locations: s_{1,2} = -zeta*omega_n +/- omega_n * sqrt(zeta^2 - 1)
 *   zeta = 0     : undamped (purely imaginary poles, sustained oscillation)
 *   0 < zeta < 1 : underdamped (complex conjugate poles, oscillatory decay)
 *   zeta = 1     : critically damped (real double pole, fastest non-oscillatory)
 *   zeta > 1     : overdamped (two real distinct poles, slow decay)
 *   zeta < 0     : unstable (poles in RHP, exponential growth)
 */
typedef struct {
    double K;       /**< DC gain */
    double zeta;    /**< damping ratio zeta */
    double omega_n; /**< natural frequency omega_n [rad/s] */
} SecondOrderModel;

/**
 * General n-th order transfer function:
 *   G(s) = (b0 + b1*s + ... + bm*s^m) / (a0 + a1*s + ... + an*s^n)
 *
 * Properness requires m <= n (den_deg >= num_deg).
 * The leading denominator coefficient a[n] should be non-zero.
 */
typedef struct {
    int     num_deg;  /**< numerator polynomial degree m */
    int     den_deg;  /**< denominator polynomial degree n (n >= m) */
    double *num;      /**< numerator coefficients [b0, b1, ..., bm] */
    double *den;      /**< denominator coefficients [a0, a1, ..., an] */
} TransferFunction;

/**
 * State-space model: dx/dt = A*x + B*u,  y = C*x + D*u
 *
 * x in R^n -- state vector
 * u in R^m -- input vector
 * y in R^p -- output vector
 *
 * Impulse response (D=0): y(t) = C * exp(A*t) * B
 * Step response (D=0):    y(t) = C * A^{-1} * (exp(A*t) - I) * B
 *
 * Matrices stored in row-major order.
 */
typedef struct {
    int     n_states;   /**< state dimension n */
    int     n_inputs;   /**< input dimension m */
    int     n_outputs;  /**< output dimension p */
    double *A;           /**< system matrix [n*n], row-major */
    double *B;           /**< input matrix [n*m], row-major */
    double *C;           /**< output matrix [p*n], row-major */
    double *D;           /**< feedthrough matrix [p*m], row-major */
} StateSpaceModel;

/* ==========================================================================
 * L1: Response Data Containers
 * ========================================================================== */

/** Single point in a response trajectory */
typedef struct {
    double t;  /**< time [s] */
    double y;  /**< output value at time t */
} ResponsePoint;

/**
 * Full response trajectory: uniformly sampled from t=0 to t=t_final.
 * N points at interval dt, where dt = t_final / (N-1).
 */
typedef struct {
    int            n_points;  /**< number of sample points */
    double         dt;        /**< uniform sampling interval [s] */
    double         t_final;   /**< final simulation time [s] */
    ResponsePoint *data;      /**< array of n_points response points */
} ResponseTrajectory;

/**
 * Time-domain performance metrics (Ogata Ch.5, Dorf Ch.4).
 *
 * All times in seconds. Overshoot stored as fraction (not percent).
 * Error integrals (IAE, ISE, ITAE) computed relative to final steady-state.
 */
typedef struct {
    double rise_time;         /**< t_r: 10%-90% rise time [s] */
    double peak_time;         /**< t_p: time to first peak [s] */
    double peak_value;        /**< y(t_p): maximum value */
    double overshoot;         /**< M_p = (y_peak - y_ss)/y_ss (fraction) */
    double overshoot_pct;     /**< overshoot in percent */
    double settling_time;     /**< t_s: time to stay within band of steady-state [s] */
    double delay_time;        /**< t_d: time to reach 50% of steady-state [s] */
    double steady_state;      /**< y_ss = lim_{t->inf} y(t) */
    double settling_band;     /**< tolerance band used (e.g., 0.02 for 2%) */
    int    n_oscillations;    /**< number of oscillations before settling */
    double decay_ratio;       /**< ratio of second peak to first peak */
    double integral_abs_error;     /**< IAE = integral |e(t)| dt */
    double integral_sq_error;      /**< ISE = integral [e(t)]^2 dt */
    double integral_time_abs_error; /**< ITAE = integral t*|e(t)| dt */
} ResponseMetrics;

/* ==========================================================================
 * L2: Common Enumerations
 * ========================================================================== */

/** Damping classification for second-order systems */
typedef enum {
    DAMPING_UNDAMPED          = 0,  /**< zeta = 0: sustained oscillation */
    DAMPING_UNDERDAMPED       = 1,  /**< 0 < zeta < 1: oscillatory decay */
    DAMPING_CRITICALLY_DAMPED = 2,  /**< zeta = 1: fastest non-oscillatory */
    DAMPING_OVERDAMPED        = 3,  /**< zeta > 1: slow non-oscillatory */
    DAMPING_UNSTABLE          = 4   /**< zeta < 0: exponential growth */
} DampingType;

/** Input signal types for time-domain simulation */
typedef enum {
    INPUT_IMPULSE,       /**< Dirac delta delta(t) */
    INPUT_STEP,          /**< Heaviside step u(t) */
    INPUT_RAMP,          /**< Ramp r(t) = t*u(t) */
    INPUT_PARABOLIC,     /**< Parabolic p(t) = (t^2/2)*u(t) */
    INPUT_SINUSOIDAL,    /**< Sinusoidal A*sin(omega*t+phi)*u(t) */
    INPUT_RECTANGULAR,   /**< Rectangular pulse of width T */
    INPUT_ARBITRARY      /**< User-defined signal sample array */
} InputSignalType;

/** Numerical quadrature method for convolution integration */
typedef enum {
    QUAD_RECTANGLE_LEFT,   /**< Left Riemann sum, O(h) */
    QUAD_RECTANGLE_RIGHT,  /**< Right Riemann sum, O(h) */
    QUAD_TRAPEZOIDAL,      /**< Trapezoidal rule, O(h^2) */
    QUAD_SIMPSON,          /**< Simpson's 1/3 rule, O(h^4) */
    QUAD_SIMPSON_38,       /**< Simpson's 3/8 rule, O(h^4) */
    QUAD_BOOLE             /**< Boole's rule (Newton-Cotes 4), O(h^6) */
} QuadMethod;

/* ==========================================================================
 * L2: Function Declarations -- System Model Utilities
 * ========================================================================== */

/**
 * Classify the damping type of a second-order model.
 * @param sys  Pointer to initialized SecondOrderModel.
 * @return DampingType enum value.
 * Complexity: O(1)
 */
DampingType classify_damping(const SecondOrderModel *sys);

/**
 * Compute the characteristic poles of a second-order system.
 * @param sys  System model.
 * @param s1   Output: first pole as [real, imag].
 * @param s2   Output: second pole as [real, imag].
 *
 * Formula: s_{1,2} = -zeta*omega_n +/- omega_n * sqrt(zeta^2 - 1)
 *   For zeta >= 1: both poles real, imag parts zero.
 *   For 0 < zeta < 1: complex conjugate pair.
 * Complexity: O(1)
 */
void second_order_poles(const SecondOrderModel *sys, double s1[2], double s2[2]);

/**
 * Compute damping ratio from percent overshoot.
 * @param overshoot_pct  Percent overshoot (e.g., 20.0 for 20%).
 * @return Damping ratio zeta.
 *
 * Formula: zeta = -ln(%OS/100) / sqrt(pi^2 + ln^2(%OS/100))
 * Textbook: Ogata Section 5-3, equation (5-18).
 * Complexity: O(1)
 */
double damping_from_overshoot(double overshoot_pct);

/**
 * Compute natural frequency from settling time and damping ratio.
 * @param ts     Settling time (2% criterion) [s].
 * @param zeta   Damping ratio.
 * @return Natural frequency omega_n [rad/s].
 *
 * Formula: omega_n = 4 / (zeta * ts)  (for 2% settling band)
 * Complexity: O(1)
 */
double omega_n_from_settling_time(double ts, double zeta);

/**
 * Compute DC gain from steady-state step response.
 * @param y_ss            Steady-state output value.
 * @param input_magnitude Magnitude of step input (default 1.0 for unit step).
 * @return DC gain K.
 * Complexity: O(1)
 */
double dc_gain_from_steady_state(double y_ss, double input_magnitude);

/**
 * Evaluate transfer function G(s) at complex frequency s.
 * @param tf   Transfer function.
 * @param s    Complex frequency [real, imag].
 * @param out  Output: G(s) as [real, imag].
 *
 * Evaluates numerator polynomial / denominator polynomial at s.
 * Returns INFINITY in both components if denominator evaluates to zero.
 * Complexity: O(num_deg + den_deg)
 */
void transfer_function_eval(const TransferFunction *tf, const double s[2], double out[2]);

/**
 * Compute impulse response of a state-space model using matrix exponential.
 * @param sys      State-space model (D=0 assumed for impulse).
 * @param t_final  End time [s].
 * @param n_steps  Number of simulation steps.
 * @return Allocated ResponseTrajectory (caller must response_trajectory_free).
 *
 * Method: y(t) = C * exp(A*t) * B
 * Matrix exponential via scaling-and-squaring + Pade(6,6).
 * Reference: Moler & Van Loan, "Nineteen Dubious Ways..." (2003), method 3.
 * Complexity: O(n^3 * log2(||A||)) per step
 */
ResponseTrajectory *state_space_impulse_response(const StateSpaceModel *sys,
                                                  double t_final, int n_steps);

/**
 * Free a ResponseTrajectory and its data array.
 */
void response_trajectory_free(ResponseTrajectory *rt);

/**
 * Create a uniformly spaced time vector from 0 to t_final.
 * @param t_final  Final time [s].
 * @param n_points Number of points (must be >= 2).
 * @return Allocated double array of n_points times (caller frees).
 * Complexity: O(n_points)
 */
double *create_time_vector(double t_final, int n_points);

/**
 * Extract a time segment from a response trajectory.
 * @param src     Source trajectory.
 * @param t_start Start time (clamped to [0, src->t_final]).
 * @param t_end   End time.
 * @return New allocated trajectory segment (caller frees).
 * Complexity: O(n_points)
 */
ResponseTrajectory *response_trajectory_slice(const ResponseTrajectory *src,
                                               double t_start, double t_end);

/**
 * Convert a second-order model to state-space (controllable canonical form).
 *   A = [[0, 1], [-omega_n^2, -2*zeta*omega_n]]
 *   B = [[0], [K * omega_n^2]]
 *   C = [[1, 0]]
 *   D = [[0]]
 *
 * @param sys  Second-order model.
 * @param ss   Output: state-space model (caller must state_space_free_matrices).
 * Complexity: O(1)
 */
void second_order_to_state_space(const SecondOrderModel *sys, StateSpaceModel *ss);

/**
 * Free all dynamically allocated matrices inside a StateSpaceModel.
 */
void state_space_free_matrices(StateSpaceModel *ss);

/**
 * Free all dynamically allocated arrays inside a TransferFunction.
 */
void transfer_function_free(TransferFunction *tf);

/**
 * Create a deep copy of a TransferFunction.
 * @return Newly allocated copy (caller must transfer_function_free).
 */
TransferFunction *transfer_function_clone(const TransferFunction *tf);

/**
 * Convert first-order model to TransferFunction.
 * G(s) = K / (tau*s + 1)
 * num = [K], den = [1, tau]
 */
TransferFunction *first_order_to_tf(const FirstOrderModel *fo);

/**
 * Convert second-order model to TransferFunction.
 * G(s) = K*omega_n^2 / (s^2 + 2*zeta*omega_n*s + omega_n^2)
 * num = [K*omega_n^2], den = [omega_n^2, 2*zeta*omega_n, 1]
 */
TransferFunction *second_order_to_tf(const SecondOrderModel *so);

/**
 * Normalize polynomial so leading coefficient is 1.0.
 * @param coeffs  Coefficient array (modified in place).
 * @param deg     Degree (index of leading coefficient).
 * @return The scaling factor applied (leading coefficient before normalization).
 * Complexity: O(deg)
 */
double polynomial_normalize(double *coeffs, int deg);

/**
 * Evaluate polynomial at real s using Horner's method.
 * P(s) = coeffs[0] + coeffs[1]*s + ... + coeffs[deg]*s^deg
 * Complexity: O(deg)
 */
double polynomial_eval_real(const double *coeffs, int deg, double s);

/**
 * Evaluate polynomial at complex s using Horner's method.
 * Complexity: O(deg)
 */
void polynomial_eval_complex(const double *coeffs, int deg,
                             const double s[2], double out[2]);

/**
 * n-th derivative of polynomial at s=0: f^{(n)}(0) = n! * coeffs[n].
 * Returns 0.0 if n > deg.
 * Complexity: O(1)
 */
double polynomial_derivative_at_zero(const double *coeffs, int deg, int n);

/**
 * Polynomial long division: num/den = quotient + remainder/den.
 * @param num       Numerator coefficients [0..num_deg].
 * @param num_deg   Numerator degree.
 * @param den       Denominator coefficients [0..den_deg].
 * @param den_deg   Denominator degree.
 * @param quotient  Output: quotient coefficients (size num_deg-den_deg+1).
 * @param remainder Output: remainder coefficients (size den_deg).
 * Complexity: O(num_deg * den_deg)
 */
void polynomial_division(const double *num, int num_deg,
                          const double *den, int den_deg,
                          double *quotient, double *remainder);

/**
 * Static gain of transfer function: K = num[0] / den[0].
 * Returns INFINITY if den[0] == 0 (system has pole at origin).
 */
double transfer_function_static_gain(const TransferFunction *tf);

/**
 * Check if transfer function is strictly proper (num_deg < den_deg).
 * Complexity: O(1)
 */
int transfer_function_is_strictly_proper(const TransferFunction *tf);

/**
 * Check BIBO stability via Routh-Hurwitz criterion.
 * @return 1 if stable (all poles in LHP), 0 if unstable or marginally stable.
 * Complexity: O(n^2) where n = den_deg
 */
int transfer_function_is_stable(const TransferFunction *tf);

/**
 * Dominant time constant = -1/Re(p_max) where p_max is rightmost pole.
 * Returns INFINITY if system is marginally stable or unstable.
 */
double dominant_time_constant(const TransferFunction *tf);

/* ==========================================================================
 * L3: Basic Mathematical Utilities
 * ========================================================================== */

/** Factorial n! (exact for n <= 20 via integer, approximate for n > 20) */
double factorial(int n);

/** Binomial coefficient C(n,k) = n!/(k!*(n-k)!) */
double binomial(int n, int k);

/** Matrix-vector multiplication: y = A*x, A is n-by-n, x and y length n */
void mat_vec_mul(const double *A, const double *x, double *y, int n);

/** Matrix-matrix multiplication: C = A*B, all n-by-n */
void mat_mat_mul(const double *A, const double *B, double *C, int n);

/**
 * Matrix exponential exp(A*t) via scaling-and-squaring with Pade(6,6).
 * @param A     n-by-n matrix (row-major).
 * @param n     Matrix dimension.
 * @param t     Time scalar.
 * @param expAt Output: n-by-n matrix exp(A*t) (row-major).
 * Reference: Higham (2005), "The Scaling and Squaring Method..."
 * Complexity: O(n^3 * log2(||A*t||))
 */
void matrix_exponential(const double *A, int n, double t, double *expAt);

/**
 * Solve Ax = b via Gaussian elimination with partial pivoting.
 * @param A  n-by-n matrix (modified in place).
 * @param b  Right-hand side vector (modified in place).
 * @param x  Output: solution vector length n.
 * @param n  Dimension.
 * @return 0 on success, -1 if singular.
 * Complexity: O(n^3)
 */
int linear_solve(const double *A, const double *b, double *x, int n);

/**
 * Eigenvalues of 2x2 real matrix [[a, b], [c, d]].
 * Uses quadratic formula on characteristic polynomial.
 */
void eigenvalues_2x2(double a, double b, double c, double d,
                     double *re1, double *im1, double *re2, double *im2);

/** LU decomposition with partial pivoting, in-place on A.
 *  Returns permutation vector piv (size n, caller-allocated). */
int lu_decompose(double *A, int n, int *piv);

/** Solve Ax = b given LU decomposition. */
void lu_solve(const double *LU, const int *piv, const double *b, double *x, int n);

/** Matrix inverse via LU decomposition. A_inv overwrites A. */
int matrix_inverse(double *A, int n);

/** Vector 2-norm. */
double vector_norm(const double *v, int n);

/** Matrix infinity norm (max row sum of absolute values). */
double matrix_inf_norm(const double *A, int n);

/** Matrix 1-norm (max column sum of absolute values). */
double matrix_one_norm(const double *A, int n);

#endif /* TIME_RESPONSE_COMMON_H */
