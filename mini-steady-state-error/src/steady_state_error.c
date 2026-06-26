/* ============================================================================
 * mini-steady-state-error: Core Steady-State Error Implementation
 * L1-L5: Full SSE analysis pipeline
 * Reference: Nise Ch.7, Ogata Ch.5, Franklin Ch.4
 * ============================================================================ */
#include "steady_state_error.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/* Memory Management */
TransferFunction tf_alloc(int numer_order, int denom_order)
{
    TransferFunction tf;
    tf.numer_order = numer_order >= 0 ? numer_order : 0;
    tf.denom_order = denom_order >= 0 ? denom_order : 0;
    tf.gain = 1.0;
    int nsize = (tf.numer_order + 1) > 0 ? (tf.numer_order + 1) : 1;
    int dsize = (tf.denom_order + 1) > 0 ? (tf.denom_order + 1) : 1;
    tf.numer_coeffs = (double *)calloc((size_t)nsize, sizeof(double));
    tf.denom_coeffs = (double *)calloc((size_t)dsize, sizeof(double));
    if (!tf.numer_coeffs || !tf.denom_coeffs) {
        free(tf.numer_coeffs); free(tf.denom_coeffs);
        tf.numer_coeffs = NULL; tf.denom_coeffs = NULL;
        tf.numer_order = -1; tf.denom_order = -1;
    }
    return tf;
}
void tf_free(TransferFunction *tf)
{
    if (!tf) return;
    free(tf->numer_coeffs); free(tf->denom_coeffs);
    tf->numer_coeffs = NULL; tf->denom_coeffs = NULL;
    tf->numer_order = 0; tf->denom_order = 0; tf->gain = 0.0;
}
TransferFunction tf_copy(const TransferFunction *tf)
{
    TransferFunction cpy = tf_alloc(tf->numer_order, tf->denom_order);
    if (cpy.numer_order < 0) return cpy;
    cpy.gain = tf->gain;
    memcpy(cpy.numer_coeffs, tf->numer_coeffs,
           (size_t)(tf->numer_order + 1) * sizeof(double));
    memcpy(cpy.denom_coeffs, tf->denom_coeffs,
           (size_t)(tf->denom_order + 1) * sizeof(double));
    return cpy;
}

void tf_print(const TransferFunction *tf, const char *label)
{
    if (!tf || !tf->numer_coeffs || !tf->denom_coeffs) {
        printf("%s: (invalid TF)\n", label ? label : "TF");
        return;
    }
    printf("%s: G(s) = %.4g * (", label ? label : "TF", tf->gain);
    int i, first_num = 1;
    for (i = 0; i <= tf->numer_order; i++) {
        if (fabs(tf->numer_coeffs[i]) > 1e-15 || (i == 0 && tf->numer_order == 0)) {
            if (!first_num) printf(" + ");
            printf("%.4g", tf->numer_coeffs[i]);
            if (i == 1) printf("*s");
            else if (i > 1) printf("*s^%d", i);
            first_num = 0;
        }
    }
    if (first_num) printf("0");
    printf(") / (");
    int first_den = 1;
    for (i = 0; i <= tf->denom_order; i++) {
        if (fabs(tf->denom_coeffs[i]) > 1e-15 || (i == 0 && tf->denom_order == 0)) {
            if (!first_den) printf(" + ");
            printf("%.4g", tf->denom_coeffs[i]);
            if (i == 1) printf("*s");
            else if (i > 1) printf("*s^%d", i);
            first_den = 0;
        }
    }
    if (first_den) printf("0");
    printf(")\n");
}

void sse_analysis_print(const SSEAnalysis *analysis)
{
    if (!analysis) return;
    printf("=== Steady-State Error Analysis ===\n");
    printf("System Type: %d\n", analysis->constants.sys_type);
    printf("Error Constants: Kp=");
    if (analysis->constants.Kp_is_inf) printf("INF");
    else printf("%.6g", analysis->constants.Kp);
    printf(", Kv=");
    if (analysis->constants.Kv_is_inf) printf("INF");
    else printf("%.6g", analysis->constants.Kv);
    printf(", Ka=");
    if (analysis->constants.Ka_is_inf) printf("INF");
    else printf("%.6g", analysis->constants.Ka);
    printf("\nSteady-State Errors:\n");
    printf("  e_ss(step) = ");
    if (isinf(analysis->error.e_ss_step)) printf("INF\n");
    else printf("%.6g\n", analysis->error.e_ss_step);
    printf("  e_ss(ramp) = ");
    if (isinf(analysis->error.e_ss_ramp)) printf("INF\n");
    else printf("%.6g\n", analysis->error.e_ss_ramp);
    printf("  e_ss(parabola) = ");
    if (isinf(analysis->error.e_ss_parabola)) printf("INF\n");
    else printf("%.6g\n", analysis->error.e_ss_parabola);
    printf("Stability: %s (FVT applicable: %s)\n",
           analysis->stability.is_stable ? "STABLE" : "UNSTABLE",
           analysis->stability.fvt_applicable ? "YES" : "NO");
    printf("DC Gain: %.6g\n", analysis->dc_gain);
    printf("Tracking Error: %.6g\n", analysis->tracking_error);
    printf("Sensitivity: dE/dK_step=%.6g dE/dK_ramp=%.6g dE/dK_para=%.6g\n",
           analysis->sensitivity[0], analysis->sensitivity[1], analysis->sensitivity[2]);
}

/* ============================================================================
 * L1: Steady-State Error from Error Constants (Nise Table 7.2)
 * e_ss_step = 1/(1+Kp), e_ss_ramp = 1/Kv, e_ss_parabola = 1/Ka
 * ============================================================================ */

SteadyStateError sse_compute_from_constants(double Kp, double Kv, double Ka)
{
    SteadyStateError err;
    double huge = 1.0e308;
    bool kp_inf = isinf(Kp) || Kp > huge;
    bool kv_inf = isinf(Kv) || Kv > huge;
    bool ka_inf = isinf(Ka) || Ka > huge;

    if (kp_inf) err.e_ss_step = 0.0;
    else if (fabs(Kp + 1.0) < 1e-15) err.e_ss_step = INFINITY;
    else err.e_ss_step = 1.0 / (1.0 + Kp);

    if (kv_inf) err.e_ss_ramp = 0.0;
    else if (fabs(Kv) < 1e-15) err.e_ss_ramp = INFINITY;
    else err.e_ss_ramp = 1.0 / Kv;

    if (ka_inf) err.e_ss_parabola = 0.0;
    else if (fabs(Ka) < 1e-15) err.e_ss_parabola = INFINITY;
    else err.e_ss_parabola = 1.0 / Ka;

    return err;
}

bool sse_is_finite_error(SystemType type, TestInputType input)
{
    switch (type) {
    case SYSTEM_TYPE_0: return (input == INPUT_STEP);
    case SYSTEM_TYPE_1: return (input == INPUT_STEP || input == INPUT_RAMP);
    case SYSTEM_TYPE_2: return (input == INPUT_STEP || input == INPUT_RAMP
                                 || input == INPUT_PARABOLA);
    default: return true;
    }
}

double sse_compute_specific(SystemType type, TestInputType input,
                             double Kp, double Kv, double Ka,
                             double magnitude)
{
    (void)type;
    SteadyStateError base = sse_compute_from_constants(Kp, Kv, Ka);
    double e_ss_unit;
    switch (input) {
    case INPUT_STEP:     e_ss_unit = base.e_ss_step;    break;
    case INPUT_RAMP:     e_ss_unit = base.e_ss_ramp;    break;
    case INPUT_PARABOLA: e_ss_unit = base.e_ss_parabola;break;
    default:             return NAN;
    }
    if (isinf(e_ss_unit)) return INFINITY;
    return e_ss_unit * magnitude;
}

/* ============================================================================
 * L2: Final Value Theorem Application
 *
 * E(s) = R(s)/(1+GH) = R(s)*S(s),  e_ss = lim_{s->0} s*E(s)
 * Step: e_ss = S(0) = 1/(1+Kp)
 * Ramp: e_ss = lim S(s)/s = 1/Kv
 * Parabola: e_ss = lim S(s)/s^2 = 1/Ka
 * ============================================================================ */

double sse_final_value_theorem(const FeedbackSystem *system, TestInputType input)
{
    if (!system) return NAN;
    const TransferFunction *G = &system->forward_path;
    double Kp, Kv, Ka;

    if (system->is_unity) {
        Kp = sse_limit_s_power_g(G, 0);
        Kv = sse_limit_s_power_g(G, 1);
        Ka = sse_limit_s_power_g(G, 2);
    } else {
        double G0 = sse_dc_gain(G);
        double H0 = sse_dc_gain(&system->feedback_path);
        Kp = G0 * H0;
        double lim_sG = sse_limit_s_power_g(G, 1);
        double lim_sH = sse_limit_s_power_g(&system->feedback_path, 1);
        Kv = lim_sG * H0 + G0 * lim_sH;
        if (fabs(Kv) > 1e308) Kv = INFINITY;
        Ka = INFINITY;
    }

    SteadyStateError err = sse_compute_from_constants(Kp, Kv, Ka);
    switch (input) {
    case INPUT_STEP:     return err.e_ss_step;
    case INPUT_RAMP:     return err.e_ss_ramp;
    case INPUT_PARABOLA: return err.e_ss_parabola;
    default:             return NAN;
    }
}

double sse_fvt_generic(const TransferFunction *error_tf)
{
    if (!error_tf || !error_tf->numer_coeffs || !error_tf->denom_coeffs)
        return NAN;

    int n_first = 0;
    while (n_first <= error_tf->numer_order &&
           fabs(error_tf->numer_coeffs[n_first]) < 1e-15) n_first++;
    int d_first = 0;
    while (d_first <= error_tf->denom_order &&
           fabs(error_tf->denom_coeffs[d_first]) < 1e-15) d_first++;

    int eff_num_first = n_first + 1;
    if (eff_num_first > d_first) return 0.0;
    else if (eff_num_first < d_first) return INFINITY;
    else {
        double nn = error_tf->numer_coeffs[n_first];
        double dd = error_tf->denom_coeffs[d_first];
        if (fabs(dd) < 1e-15) return NAN;
        return error_tf->gain * nn / dd;
    }
}

/* ============================================================================
 * L3: Transfer Function Evaluation at Limits
 * ============================================================================ */

double sse_dc_gain(const TransferFunction *tf)
{
    if (!tf || !tf->numer_coeffs || !tf->denom_coeffs) return NAN;
    if (tf->numer_order < 0 || tf->denom_order < 0) return NAN;
    double n0 = tf->numer_coeffs[0];
    double d0 = tf->denom_coeffs[0];
    if (fabs(d0) < 1e-15) {
        if (fabs(n0) < 1e-15) return NAN;
        return (n0 > 0) ? INFINITY : -INFINITY;
    }
    return tf->gain * n0 / d0;
}

double sse_limit_s_power_g(const TransferFunction *tf, int k)
{
    if (!tf || !tf->numer_coeffs || !tf->denom_coeffs) return NAN;
    if (k < 0) return NAN;

    /* Find first non-zero coefficient in numerator and denominator.
     * lim_{s->0} s^k * N(s)/D(s):
     * Let n_first = min index with non-zero numer coeff
     * Let d_first = min index with non-zero denom coeff
     * Effective numerator start after *s^k: k + n_first
     * If k+n_first == d_first: limit = n_{n_first}/d_{d_first}
     * If k+n_first > d_first:  limit = 0
     * If k+n_first < d_first:  limit = infinity */

    int n_first = 0;
    while (n_first <= tf->numer_order &&
           fabs(tf->numer_coeffs[n_first]) < 1e-15) n_first++;
    if (n_first > tf->numer_order) n_first = tf->numer_order + 1;

    int d_first = 0;
    while (d_first <= tf->denom_order &&
           fabs(tf->denom_coeffs[d_first]) < 1e-15) d_first++;
    if (d_first > tf->denom_order) return NAN;

    int eff_num_first = k + n_first;

    if (eff_num_first > d_first) return 0.0;
    else if (eff_num_first < d_first) return INFINITY;
    else {
        double n_val = tf->numer_coeffs[n_first];
        double d_val = tf->denom_coeffs[d_first];
        if (fabs(d_val) < 1e-15) return NAN;
        return tf->gain * n_val / d_val;
    }
}

/* ============================================================================
 * L2: Error Constants and System Type from Transfer Function
 * ============================================================================ */

ErrorConstants sse_error_constants_from_tf(const TransferFunction *open_loop_G)
{
    ErrorConstants ec;
    memset(&ec, 0, sizeof(ec));
    if (!open_loop_G) return ec;

    SystemType st = sse_determine_system_type(open_loop_G);
    ec.sys_type = st;

    double kp = sse_limit_s_power_g(open_loop_G, 0);
    double kv = sse_limit_s_power_g(open_loop_G, 1);
    double ka = sse_limit_s_power_g(open_loop_G, 2);

    ec.Kp_is_inf = isinf(kp);
    ec.Kv_is_inf = isinf(kv);
    ec.Ka_is_inf = isinf(ka);
    ec.Kp = ec.Kp_is_inf ? 0.0 : kp;
    ec.Kv = ec.Kv_is_inf ? 0.0 : kv;
    ec.Ka = ec.Ka_is_inf ? 0.0 : ka;

    return ec;
}

SystemType sse_determine_system_type(const TransferFunction *tf)
{
    if (!tf || !tf->denom_coeffs) return SYSTEM_TYPE_0;
    int type = 0;
    while (type <= tf->denom_order &&
           fabs(tf->denom_coeffs[type]) < 1e-15) type++;
    if (type == 0) return SYSTEM_TYPE_0;
    if (type == 1) return SYSTEM_TYPE_1;
    if (type == 2) return SYSTEM_TYPE_2;
    return SYSTEM_TYPE_N;
}

/* ============================================================================
 * L4: Routh-Hurwitz Stability Criterion (Ogata 5.6, Nise 6.2)
 *
 * Routh array from polynomial a0 + a1*s + ... + a_n*s^n.
 * Row 0 (s^n):   a_n, a_{n-2}, a_{n-4}, ...
 * Row 1 (s^{n-1}): a_{n-1}, a_{n-3}, a_{n-5}, ...
 * Row r (r>=2): b_j = (p*q - r*s)/p where
 *   p = routh[r-2][0], q = routh[r-2][j+1],
 *   r = routh[r-1][0], s = routh[r-1][j+1]
 *
 * Stability condition: all elements in first column have same sign.
 * Number of RHP poles = number of sign changes in first column.
 *
 * Special cases handled:
 * 1. Zero in first column -> epsilon method (replace zero with small epsilon)
 * 2. Entire row zero -> use coefficients of auxiliary polynomial derivative
 * ============================================================================ */

StabilityCheck sse_routh_hurwitz_check(const double *closed_loop_denom, int order)
{
    StabilityCheck result;
    memset(&result, 0, sizeof(result));

    if (!closed_loop_denom || order <= 0) {
        result.unstable_pole_count = -1;
        return result;
    }

    int cols = (order + 2) / 2;
    double **routh = (double **)malloc((size_t)(order + 1) * sizeof(double *));
    int r, c;
    for (r = 0; r <= order; r++) {
        routh[r] = (double *)calloc((size_t)cols, sizeof(double));
    }

    /* Row 0: even powers (a_n, a_{n-2}, a_{n-4}, ...) */
    int col0 = 0;
    for (c = order; c >= 0; c -= 2) routh[0][col0++] = closed_loop_denom[c];

    /* Row 1: odd powers (a_{n-1}, a_{n-3}, a_{n-5}, ...) */
    int col1 = 0;
    for (c = order - 1; c >= 0; c -= 2) routh[1][col1++] = closed_loop_denom[c];

    int sign_changes = 0;
    double prev_sign = 0.0;
    bool prev_sign_set = false;
    double epsilon = 1e-9;

    for (r = 2; r <= order; r++) {
        int max_col = (order - r + 2) / 2;
        if (max_col < 1) max_col = 1;

        for (c = 0; c < max_col - 1; c++) {
            double a = routh[r-2][0];   /* first element of row r-2 */
            double b = routh[r-2][c+1]; /* element at column c+1 of row r-2 */
            double d = routh[r-1][0];   /* first element of row r-1 */
            double e = routh[r-1][c+1]; /* element at column c+1 of row r-1 */

            if (fabs(d) < 1e-15) d = epsilon; /* epsilon method */
            routh[r][c] = (a * e - d * b) / d;
        }

        /* Check if entire row is zero -> auxiliary polynomial case */
        bool entire_row_zero = true;
        for (c = 0; c < cols; c++) {
            if (fabs(routh[r][c]) > 1e-12) {
                entire_row_zero = false;
                break;
            }
        }

        if (entire_row_zero && r < order) {
            /* Entire row zero: use derivative of auxiliary polynomial.
             * Auxiliary polynomial is formed from the row above (r-1).
             * Powers: s^{order-r+2}, s^{order-r}, s^{order-r-2}, ...
             * Derivative multiplies each term by its power. */
            int power = order - r + 2;
            for (c = 0; c < cols && power > 0; c++) {
                routh[r][c] = power * routh[r-1][c];
                power -= 2;
            }
            while (c < cols) { routh[r][c] = 0.0; c++; }
        }
    }

    /* Count sign changes in the first column */
    for (r = 0; r <= order; r++) {
        double val = routh[r][0];

        /* Handle zero or near-zero entries with epsilon method */
        if (fabs(val) < 1e-12) {
            val = epsilon;
            epsilon = -epsilon; /* flip sign to avoid bias */
        }

        if (!prev_sign_set) {
            prev_sign = val;
            prev_sign_set = true;
        } else {
            if ((prev_sign > 0 && val < -1e-12) ||
                (prev_sign < -1e-12 && val > 0)) {
                sign_changes++;
            }
            if (fabs(val) > 1e-12) prev_sign = val;
        }
    }

    result.is_stable = (sign_changes == 0);
    result.unstable_pole_count = sign_changes;
    result.fvt_applicable = result.is_stable;
    result.dominant_pole_real = 0.0;

    for (r = 0; r <= order; r++) free(routh[r]);
    free(routh);
    return result;
}

/* ============================================================================
 * L5: Full Steady-State Error Analysis Pipeline
 *
 * Computes error constants, steady-state errors for all standard inputs,
 * checks stability via Routh-Hurwitz, and computes gain sensitivity.
 * ============================================================================ */

SSEAnalysis sse_full_analysis(const FeedbackSystem *system)
{
    SSEAnalysis analysis;
    memset(&analysis, 0, sizeof(analysis));
    if (!system) return analysis;

    /* Step 1: Compute open-loop G_ol = G*H */
    TransferFunction open_loop;
    const TransferFunction *G = &system->forward_path;
    const TransferFunction *H = &system->feedback_path;

    if (system->is_unity) {
        open_loop = tf_copy(G);
    } else {
        int no = G->numer_order + H->numer_order;
        int do_ = G->denom_order + H->denom_order;
        open_loop = tf_alloc(no, do_);
        if (open_loop.numer_order >= 0) {
            open_loop.gain = G->gain * H->gain;
            int i, j;
            for (i = 0; i <= G->numer_order; i++)
                for (j = 0; j <= H->numer_order; j++)
                    open_loop.numer_coeffs[i+j] +=
                        G->numer_coeffs[i] * H->numer_coeffs[j];
            for (i = 0; i <= G->denom_order; i++)
                for (j = 0; j <= H->denom_order; j++)
                    open_loop.denom_coeffs[i+j] +=
                        G->denom_coeffs[i] * H->denom_coeffs[j];
        }
    }

    /* Step 2: Compute error constants */
    analysis.constants = sse_error_constants_from_tf(&open_loop);

    /* Step 3: Compute steady-state errors for standard inputs */
    analysis.error = sse_compute_from_constants(
        analysis.constants.Kp_is_inf ? INFINITY : analysis.constants.Kp,
        analysis.constants.Kv_is_inf ? INFINITY : analysis.constants.Kv,
        analysis.constants.Ka_is_inf ? INFINITY : analysis.constants.Ka);

    /* Step 4: Compute closed-loop transfer function */
    TransferFunction cl = sse_closed_loop_tf(system);

    /* Step 5: Stability check via Routh-Hurwitz on closed-loop denominator */
    if (cl.denom_order >= 0 && cl.denom_coeffs) {
        analysis.stability = sse_routh_hurwitz_check(
            cl.denom_coeffs, cl.denom_order);
    }

    /* Step 6: DC gain and tracking error */
    analysis.dc_gain = sse_dc_gain(&cl);
    analysis.tracking_error = 1.0 - analysis.dc_gain;

    /* Step 7: Gain sensitivity analysis
     * For Type 0 + step:  de_ss/dK = -1/(1+K)^2 = -e_ss^2
     * For Type 1 + ramp:  de_ss/dKv = -1/Kv^2 = -e_ss/Kv
     * For Type 2 + parab: de_ss/dKa = -1/Ka^2 = -e_ss/Ka */
    if (!isinf(analysis.error.e_ss_step) && analysis.error.e_ss_step > 0)
        analysis.sensitivity[0] = -analysis.error.e_ss_step
                                   * analysis.error.e_ss_step;
    if (!isinf(analysis.error.e_ss_ramp) && analysis.error.e_ss_ramp > 0
        && !analysis.constants.Kv_is_inf)
        analysis.sensitivity[1] = -analysis.error.e_ss_ramp
                                   / analysis.constants.Kv;
    if (!isinf(analysis.error.e_ss_parabola) && analysis.error.e_ss_parabola > 0
        && !analysis.constants.Ka_is_inf)
        analysis.sensitivity[2] = -analysis.error.e_ss_parabola
                                   / analysis.constants.Ka;

    tf_free(&open_loop);
    tf_free(&cl);
    return analysis;
}

/* ============================================================================
 * L1: Closed-Loop Transfer Function Computation
 *
 * T(s) = G(s) / (1 + G(s)H(s))
 *
 * Let G(s) = N_G/D_G, H(s) = N_H/D_H.
 * T(s) = (N_G * D_H) / (D_G * D_H + N_G * N_H)
 *
 * For unity feedback (H=1): T(s) = N_G / (D_G + N_G)
 * ============================================================================ */

TransferFunction sse_closed_loop_tf(const FeedbackSystem *system)
{
    TransferFunction result = {NULL, -1, NULL, -1, 0.0};
    if (!system) return result;

    const TransferFunction *G = &system->forward_path;
    const TransferFunction *H = &system->feedback_path;

    if (system->is_unity) {
        /* T(s) = N_G(s) / (D_G(s) + N_G(s)) */
        int max_order = (G->numer_order > G->denom_order)
                         ? G->numer_order : G->denom_order;
        result = tf_alloc(G->numer_order, max_order);
        if (result.numer_order < 0) return result;
        result.gain = 1.0;

        int i;
        /* Numerator: N_G (with gain) */
        for (i = 0; i <= G->numer_order; i++)
            result.numer_coeffs[i] = G->numer_coeffs[i] * G->gain;
        /* Denominator: D_G + N_G (with gain) */
        for (i = 0; i <= G->denom_order; i++)
            result.denom_coeffs[i] = G->denom_coeffs[i];
        for (i = 0; i <= G->numer_order; i++)
            result.denom_coeffs[i] += G->gain * G->numer_coeffs[i];
    } else {
        /* General case: T(s) = (N_G*D_H) / (D_G*D_H + N_G*N_H) */
        int num_order = G->numer_order + H->denom_order;
        int den_order1 = G->denom_order + H->denom_order;
        int den_order2 = G->numer_order + H->numer_order;
        int den_order = (den_order1 > den_order2) ? den_order1 : den_order2;

        result = tf_alloc(num_order, den_order);
        if (result.numer_order < 0) return result;
        result.gain = 1.0;

        int i, j;
        /* Numerator: N_G * D_H */
        for (i = 0; i <= G->numer_order; i++)
            for (j = 0; j <= H->denom_order; j++)
                result.numer_coeffs[i+j] += G->gain
                    * G->numer_coeffs[i] * H->denom_coeffs[j];

        /* Denominator: D_G*D_H + G_gain*N_G*N_H */
        for (i = 0; i <= G->denom_order; i++)
            for (j = 0; j <= H->denom_order; j++)
                result.denom_coeffs[i+j] +=
                    G->denom_coeffs[i] * H->denom_coeffs[j];
        for (i = 0; i <= G->numer_order; i++)
            for (j = 0; j <= H->numer_order; j++)
                result.denom_coeffs[i+j] += G->gain
                    * G->numer_coeffs[i] * H->numer_coeffs[j];
    }

    return result;
}

/* ============================================================================
 * L4: Simplified FVT Precondition Check
 *
 * The Final Value Theorem requires that s*F(s) has all poles in the
 * left half-plane, except possibly a simple pole at s=0.
 *
 * For a polynomial denominator D(s), we check:
 * 1. All poles must have Re(s) <= 0
 * 2. Any pole on the imaginary axis must be simple (multiplicity 1)
 * 3. At most one simple pole at s=0
 *
 * This implementation uses the Routh-Hurwitz criterion on a shifted
 * polynomial to detect RHP poles. For imaginary-axis poles, we check
 * the Routh array for special patterns (entire row zero).
 *
 * For a simplified but reliable check, we use the Routh-Hurwitz test
 * combined with a check for zero constant term (indicating an s=0 factor).
 * ============================================================================ */

bool sse_check_fvt_precondition(const TransferFunction *error_tf, double tol)
{
    if (!error_tf || !error_tf->denom_coeffs) return false;
    int n = error_tf->denom_order;
    if (n <= 0) return true;

    /* Count the multiplicity of s=0 factor in denominator */
    int s_power = 0;
    while (s_power <= n &&
           fabs(error_tf->denom_coeffs[s_power]) < tol) {
        s_power++;
    }

    /* FVT allows at most one pole at s=0 */
    if (s_power > 1) return false;

    /* If there is an s=0 factor, remove it for stability check */
    if (s_power > 0) {
        /* Reduced denominator: D(s) = s * D'(s)
         * D'(s) coefficients are d_1, d_2, ..., d_n */
        int reduced_order = n - s_power;
        if (reduced_order <= 0) return true;

        double *reduced = (double *)malloc((size_t)(reduced_order + 1) * sizeof(double));
        int i;
        for (i = 0; i <= reduced_order; i++) {
            reduced[i] = error_tf->denom_coeffs[i + s_power];
        }

        StabilityCheck sc = sse_routh_hurwitz_check(reduced, reduced_order);
        free(reduced);
        return sc.is_stable;
    }

    /* No s=0 factor: check full denominator stability */
    StabilityCheck sc = sse_routh_hurwitz_check(
        error_tf->denom_coeffs, n);
    return sc.is_stable;
}
