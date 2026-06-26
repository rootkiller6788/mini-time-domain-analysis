#ifndef STEADY_STATE_ERROR_H
#define STEADY_STATE_ERROR_H

#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* ============================================================================
 * mini-steady-state-error: Steady-State Error Analysis for Control Systems
 *
 * Reference: Norman S. Nise, "Control Systems Engineering", 8th Ed, Ch.7
 * Reference: Ogata, "Modern Control Engineering", 5th Ed, Ch.5
 * Reference: Franklin, Powell, Emami-Naeini, "Feedback Control", 8th Ed, Ch.4
 *
 * L1 - Core Definitions: steady-state error, error constants, system type
 * L2 - Core Concepts:   Final Value Theorem, unity feedback error analysis
 * L3 - Math Structures: Laplace-domain transfer functions, polynomial ratios
 * L4 - Fundamental Laws: Final Value Theorem, stability precondition
 * L5 - Engineering Methods: Routh-Hurwitz, error constant computation, sensitivity
 * ============================================================================ */

/* L1: Steady-State Error Definition (Nise 7.2, Ogata 5.1)
 *   e_ss = lim_{t→∞} e(t) = lim_{s→0} s*E(s)   [Final Value Theorem]
 * For unity feedback G(s):
 *   e_ss_step  = 1 / (1 + Kp)      position error constant Kp
 *   e_ss_ramp  = 1 / Kv            velocity error constant Kv
 *   e_ss_parab = 1 / Ka            acceleration error constant Ka
 * System type N: number of pure integrators (1/s factors) in forward path. */

typedef struct {
    double e_ss_step;
    double e_ss_ramp;
    double e_ss_parabola;
} SteadyStateError;

typedef enum {
    SYSTEM_TYPE_0 = 0,
    SYSTEM_TYPE_1 = 1,
    SYSTEM_TYPE_2 = 2,
    SYSTEM_TYPE_N = 99
} SystemType;

typedef enum {
    INPUT_STEP = 0,
    INPUT_RAMP = 1,
    INPUT_PARABOLA = 2,
    INPUT_CUSTOM = 99
} TestInputType;

/* L3: Transfer function representation */
typedef struct {
    double *numer_coeffs;
    int    numer_order;
    double *denom_coeffs;
    int    denom_order;
    double gain;
} TransferFunction;

typedef struct {
    TransferFunction forward_path;
    TransferFunction feedback_path;
    bool             is_unity;
} FeedbackSystem;

/* L1/L2: Error constant results */
typedef struct {
    double Kp;
    double Kv;
    double Ka;
    bool   Kp_is_inf;
    bool   Kv_is_inf;
    bool   Ka_is_inf;
    SystemType sys_type;
} ErrorConstants;

/* L4: Stability check result */
typedef struct {
    bool   is_stable;
    bool   fvt_applicable;
    int    unstable_pole_count;
    double dominant_pole_real;
} StabilityCheck;

/* L5: Complete SSE analysis result */
typedef struct {
    SteadyStateError  error;
    ErrorConstants    constants;
    StabilityCheck    stability;
    double            sensitivity[3];
    double            dc_gain;
    double            tracking_error;
} SSEAnalysis;

/* ---- Core API (L1-L5) ---- */

SteadyStateError sse_compute_from_constants(double Kp, double Kv, double Ka);
bool sse_is_finite_error(SystemType type, TestInputType input);
double sse_compute_specific(SystemType type, TestInputType input,
                             double Kp, double Kv, double Ka, double magnitude);
double sse_final_value_theorem(const FeedbackSystem *system, TestInputType input);
double sse_fvt_generic(const TransferFunction *error_tf);
double sse_dc_gain(const TransferFunction *tf);
double sse_limit_s_power_g(const TransferFunction *tf, int k);
ErrorConstants sse_error_constants_from_tf(const TransferFunction *open_loop_G);
SystemType sse_determine_system_type(const TransferFunction *tf);
StabilityCheck sse_routh_hurwitz_check(const double *closed_loop_denom, int order);
bool sse_check_fvt_precondition(const TransferFunction *error_tf, double tol);
SSEAnalysis sse_full_analysis(const FeedbackSystem *system);
TransferFunction sse_closed_loop_tf(const FeedbackSystem *system);

/* Memory management */
TransferFunction tf_alloc(int numer_order, int denom_order);
void tf_free(TransferFunction *tf);
TransferFunction tf_copy(const TransferFunction *tf);
void tf_print(const TransferFunction *tf, const char *label);
void sse_analysis_print(const SSEAnalysis *analysis);

#endif /* STEADY_STATE_ERROR_H */
