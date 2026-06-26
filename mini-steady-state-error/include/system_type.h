#ifndef SYSTEM_TYPE_H
#define SYSTEM_TYPE_H

#include "steady_state_error.h"

/* ============================================================================
 * mini-steady-state-error: System Type Analysis Module
 *
 * L1: System type definition (Nise 7.3)
 * L2: Relationship between open-loop pole at origin count and error behavior
 * L3: Polynomial analysis for type detection
 * L4: Internal model principle (Francis & Wonham, 1976)
 * L5: Type augmentation design methods
 *
 * System type N = number of pure integrators (1/s) in the forward path.
 * This determines which polynomial inputs can be tracked with zero error.
 * ============================================================================ */

/* L1: Extended system type information */
typedef struct {
    int    type;                    /* 0, 1, 2, or N */
    int    integrator_count;        /* number of poles at s=0 in open-loop */
    double integrator_gain;         /* product of all integrator gains */
    bool   has_differentiator;      /* true if numerator has s factor (zero at origin) */
    int    relative_degree;         /* denom_order - numer_order */
    bool   is_proper;               /* denom_order >= numer_order */
    bool   is_strictly_proper;      /* denom_order > numer_order */
} SystemTypeInfo;

/* L2: Tracking capability by system type */
typedef enum {
    TRACK_DIVERGES = 0,    /* error grows unbounded (∞) */
    TRACK_CONSTANT_ERROR,  /* finite non-zero error */
    TRACK_ZERO_ERROR       /* perfect asymptotic tracking */
} TrackingCapability;

/* ---- L1/L2: Core system type operations ---- */

/**
 * L1: Get detailed system type information from transfer function.
 * Counts poles at origin, checks for zeros at origin, computes relative degree.
 *
 * Complexity: O(n) where n = max(numer_order, denom_order)
 */
SystemTypeInfo st_analyze(const TransferFunction *tf);

/**
 * L1: Get human-readable system type name.
 * Complexity: O(1)
 */
const char *st_type_name(int type);

/**
 * L2: Determine tracking capability for a given (system_type, input_type) pair.
 * Implements the standard table (Nise Table 7.2):
 *   Type 0 + Step  → TRACK_CONSTANT_ERROR
 *   Type 0 + Ramp  → TRACK_DIVERGES
 *   Type 1 + Step  → TRACK_ZERO_ERROR
 *   Type 1 + Ramp  → TRACK_CONSTANT_ERROR
 *   Type 1 + Para  → TRACK_DIVERGES
 *   Type 2 + Step  → TRACK_ZERO_ERROR
 *   Type 2 + Ramp  → TRACK_ZERO_ERROR
 *   Type 2 + Para  → TRACK_CONSTANT_ERROR
 *
 * Complexity: O(1)
 */
TrackingCapability st_tracking_capability(int system_type, TestInputType input);

/**
 * L2: Get the maximum polynomial degree that can be tracked with zero error.
 * For Type N: can track polynomial inputs up to degree (N-1) with zero error.
 *   Type 0: degree = -1 (cannot track any with zero error)
 *   Type 1: degree = 0 (step)
 *   Type 2: degree = 1 (ramp)
 *   Type N: degree = N-1
 *
 * Complexity: O(1)
 */
int st_max_trackable_degree(int system_type);

/**
 * L2: Get the polynomial degree that produces finite non-zero error.
 * For Type N: input degree N produces finite error.
 *
 * Complexity: O(1)
 */
int st_finite_error_degree(int system_type);

/* ---- L3: Polynomial decomposition for type analysis ---- */

/**
 * L3: Count trailing zeros in polynomial coefficient array.
 * Trailing zeros = zeros at highest powers (coefficient of highest s powers).
 * Leading zeros = zeros at lowest powers (coefficient of lowest s powers).
 *
 * For denominator a0 + a1*s + a2*s^2 + ... + an*s^n:
 *   system type = count of leading zeros (a0=0 → type≥1, a0=a1=0 → type≥2, etc.)
 *
 * Complexity: O(n)
 */
int st_count_leading_zeros(const double *coeffs, int order);

/**
 * L3: Count leading non-zero coefficient index.
 * Returns smallest index i such that coeffs[i] != 0, or -1 if all zero.
 *
 * Complexity: O(n)
 */
int st_first_nonzero_index(const double *coeffs, int order);

/**
 * L3: Factor out s^N from denominator polynomial.
 * Given denom(s) = s^N * (a_N + a_{N+1}*s + ... + a_m*s^{m-N}),
 * returns the reduced polynomial a_N + a_{N+1}*s + ...
 *
 * @param denom_coeffs  input denominator coefficients
 * @param denom_order   input denominator order
 * @param N             number of s factors to remove
 * @param reduced_order output: order of reduced polynomial
 * @param reduced_coeffs output: reduced polynomial coefficients (caller frees)
 *
 * Complexity: O(m)
 */
void st_factor_out_s_power(const double *denom_coeffs, int denom_order, int N,
                            double **reduced_coeffs, int *reduced_order);

/* ---- L4: Internal Model Principle ---- */

/**
 * L4: Internal Model Principle (Francis & Wonham, 1976, Automatica 12:457-465).
 *
 * IMP states: For perfect asymptotic tracking of a reference signal, the
 * loop must contain a model of the reference signal generator (i.e., poles
 * at the reference signal's natural frequencies on the imaginary axis).
 *
 * For step tracking: need integrator (1/s) in loop.
 * For ramp tracking: need double integrator (1/s^2).
 * For sinusoidal tracking (frequency ω0): need s^2 + ω0^2 in denominator.
 *
 * This function verifies whether the open-loop TF satisfies IMP for a given
 * reference signal type.
 *
 * Complexity: O(n) where n = denom_order
 */
bool st_check_imp_satisfied(const TransferFunction *open_loop, TestInputType ref_type);

/**
 * L4: Check IMP for sinusoidal reference of given frequency.
 * The loop must contain a factor (s^2 + omega^2) in the denominator.
 *
 * Complexity: O(n)
 */
bool st_check_imp_sinusoidal(const TransferFunction *open_loop, double omega);

/**
 * L4: List required denominator factors for IMP based on reference type.
 * Returns a null-terminated array of strings describing required poles.
 * Caller must free each string and the array.
 *
 * Complexity: O(1)
 */
char **st_imp_required_factors(TestInputType ref_type);

/* ---- L5: Type augmentation ---- */

/**
 * L5: Design a lag compensator to increase system type by 1.
 * Lag compensator: G_c(s) = (s + z_c) / (s + p_c) with |z_c| > |p_c|
 *
 * Adding a lag compensator increases system type from N to N+1
 * (adds a pole near the origin, zero placed to maintain stability).
 *
 * @param original_type  current system type
 * @param z_c            compensator zero location (negative real)
 * @param p_c            compensator pole location (negative real, |p_c| << |z_c|)
 * @return               new system type after compensation
 *
 * Complexity: O(1)
 */
int st_augment_type_lag(int original_type, double z_c, double p_c);

/**
 * L5: Design PI controller parameters to achieve Type 1 from Type 0.
 * PI controller: G_c(s) = K_p + K_i/s = K_p*(s + K_i/K_p)/s
 * This adds one integrator (Type → Type+1).
 *
 * Given desired Kv (velocity error constant) and stability margin (phase margin),
 * compute K_p and K_i.
 *
 * @param desired_Kv       desired velocity error constant
 * @param plant_dc_gain    DC gain of the plant G(0)
 * @param K_p_out          output: proportional gain
 * @param K_i_out          output: integral gain
 *
 * Complexity: O(1)
 */
void st_design_pi_for_type1(double desired_Kv, double plant_dc_gain,
                             double *K_p_out, double *K_i_out);

/**
 * L5: Design PID controller parameters to achieve Type 2 from Type 0.
 * PID: G_c(s) = K_p + K_i/s + K_d*s = (K_d*s^2 + K_p*s + K_i)/s
 *
 * @param desired_Ka       desired acceleration error constant
 * @param plant_dc_gain    DC gain of the plant
 * @param K_p_out          output: proportional gain
 * @param K_i_out          output: integral gain
 * @param K_d_out          output: derivative gain
 *
 * Complexity: O(1)
 */
void st_design_pid_for_type2(double desired_Ka, double plant_dc_gain,
                              double *K_p_out, double *K_i_out, double *K_d_out);

/**
 * L2: Compute the effective system type of a compensated system.
 * When compensator G_c(s) is cascaded with plant G_p(s),
 * the combined type = type(G_c) + type(G_p).
 *
 * Complexity: O(n+m) polynomial multiplication
 */
SystemTypeInfo st_combined_type(const TransferFunction *compensator,
                                 const TransferFunction *plant);

/* ---- Utility ---- */

void st_info_print(const SystemTypeInfo *info);
const char *st_tracking_name(TrackingCapability cap);

#endif /* SYSTEM_TYPE_H */
