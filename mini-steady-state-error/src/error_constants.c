/* ============================================================================
 * mini-steady-state-error: Error Constants Implementation
 *
 * L1: Kp, Kv, Ka definitions from transfer function limits
 * L2: System type to error constant relationship
 * L3: Polynomial evaluation, Bode plot mapping
 * L5: Non-unity feedback, generalized SSE, bandwidth estimation
 *
 * Reference: Nise 7.3-7.5, Ogata 5.2-5.3
 * ============================================================================ */

#include "error_constants.h"
#include "system_type.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * L1: Direct Error Constant Computation from Transfer Function
 * ============================================================================ */

double ec_compute_Kp(const TransferFunction *open_loop)
{
    if (!open_loop) return NAN;

    /* Kp = lim_{s->0} G(s)H(s)
     * For polynomial form: G(s) = K * (n0 + n1*s + ...) / (d0 + d1*s + ...)
     * Kp = K * n0 / d0 if d0 != 0, else Kp = infinity */

    double n0 = open_loop->numer_coeffs[0];
    double d0 = open_loop->denom_coeffs[0];

    if (fabs(d0) < 1e-15) {
        if (fabs(n0) < 1e-15) return NAN;
        return INFINITY;
    }
    return open_loop->gain * n0 / d0;
}

double ec_compute_Kv(const TransferFunction *open_loop)
{
    if (!open_loop) return NAN;

    /* Kv = lim_{s->0} s * G(s)H(s)
     * If d0 != 0: Kv = 0 (since s factor in numerator gives s*G(0) -> 0)
     * If d0 = 0, d1 != 0: Kv = K * n0 / d1
     * If d0 = d1 = 0: Kv = infinity */

    double n0 = open_loop->numer_coeffs[0];
    double d0 = open_loop->denom_coeffs[0];
    double d1 = (open_loop->denom_order >= 1) ? open_loop->denom_coeffs[1] : 0.0;

    if (fabs(d0) > 1e-15) {
        /* Type 0: denominator has constant term -> Kv = 0 */
        return 0.0;
    }

    /* d0 = 0: system has at least one integrator */
    if (fabs(d1) < 1e-15) {
        /* Type >= 2: d0 = d1 = 0 -> Kv = infinity */
        if (fabs(n0) < 1e-15) return NAN;
        return INFINITY;
    }

    /* Type 1: Kv = K * n0 / d1 */
    return open_loop->gain * n0 / d1;
}

double ec_compute_Ka(const TransferFunction *open_loop)
{
    if (!open_loop) return NAN;

    /* Ka = lim_{s->0} s^2 * G(s)H(s)
     * If d0 != 0: Ka = 0
     * If d0 = 0, d1 != 0: Ka = 0 (Type 1)
     * If d0 = d1 = 0, d2 != 0: Ka = K * n0 / d2 (Type 2)
     * If d0 = d1 = d2 = 0: Ka = infinity (Type >= 3) */

    double n0 = open_loop->numer_coeffs[0];
    double d0 = open_loop->denom_coeffs[0];
    double d1 = (open_loop->denom_order >= 1) ? open_loop->denom_coeffs[1] : 0.0;
    double d2 = (open_loop->denom_order >= 2) ? open_loop->denom_coeffs[2] : 0.0;

    if (fabs(d0) > 1e-15) return 0.0;       /* Type 0 */
    if (fabs(d1) > 1e-15) return 0.0;       /* Type 1 */
    if (fabs(d2) < 1e-15) {                  /* Type >= 3 */
        if (fabs(n0) < 1e-15) return NAN;
        return INFINITY;
    }

    /* Type 2: Ka = K * n0 / d2 */
    return open_loop->gain * n0 / d2;
}

/* ============================================================================
 * L2: System Type to Error Constant Relationship
 * ============================================================================ */

void ec_finiteness_by_type(int system_type, bool finite_out[3])
{
    /* finite_out[0] = Kp finite? [1] = Kv finite? [2] = Ka finite? */
    switch (system_type) {
    case 0:
        finite_out[0] = true;
        finite_out[1] = false;
        finite_out[2] = false;
        break;
    case 1:
        finite_out[0] = false; /* Kp = inf */
        finite_out[1] = true;
        finite_out[2] = false;
        break;
    case 2:
        finite_out[0] = false;
        finite_out[1] = false;
        finite_out[2] = true;
        break;
    default: /* Type >= 3 */
        finite_out[0] = false;
        finite_out[1] = false;
        finite_out[2] = false;
        break;
    }
}

int ec_error_formula_id(SystemType type, TestInputType input)
{
    /* Returns formula identifier:
     * 0: e_ss = 1/(1+Kp)    1: e_ss = 1/Kv    2: e_ss = 1/Ka
     * 3: e_ss = 0           4: e_ss = infinity */
    switch (type) {
    case SYSTEM_TYPE_0:
        if (input == INPUT_STEP)      return 0; /* 1/(1+Kp) */
        if (input == INPUT_RAMP)      return 4; /* inf */
        if (input == INPUT_PARABOLA)  return 4; /* inf */
        break;
    case SYSTEM_TYPE_1:
        if (input == INPUT_STEP)      return 3; /* 0 */
        if (input == INPUT_RAMP)      return 1; /* 1/Kv */
        if (input == INPUT_PARABOLA)  return 4; /* inf */
        break;
    case SYSTEM_TYPE_2:
        if (input == INPUT_STEP)      return 3; /* 0 */
        if (input == INPUT_RAMP)      return 3; /* 0 */
        if (input == INPUT_PARABOLA)  return 2; /* 1/Ka */
        break;
    default:
        if (input <= (TestInputType)(type - 1)) return 3; /* 0 */
        if (input == (TestInputType)type) {
            if (type == 0) return 0;
            if (type == 1) return 1;
            if (type == 2) return 2;
            return 1; /* generic 1/K_N */
        }
        return 4; /* inf */
    }
    return -1;
}

/* ============================================================================
 * L3: Numerical Evaluation of s^k * G(s) at Limit s->0
 * ============================================================================ */

double ec_poly_at_zero(const double *coeffs, int order)
{
    if (!coeffs || order < 0) return NAN;
    return coeffs[0];
}

double ec_evaluate_limit_s_power(const double *numer, int numer_order,
                                  const double *denom, int denom_order,
                                  int k, double *finite_flag)
{
    if (!numer || !denom || !finite_flag) return NAN;

    int n_first = 0;
    while (n_first <= numer_order && fabs(numer[n_first]) < 1e-15) n_first++;
    if (n_first > numer_order) n_first = numer_order + 1;

    int d_first = 0;
    while (d_first <= denom_order && fabs(denom[d_first]) < 1e-15) d_first++;
    if (d_first > denom_order) { *finite_flag = 0.0; return NAN; }

    int eff_num_first = k + n_first;

    if (eff_num_first > d_first) {
        *finite_flag = 1.0;
        return 0.0;
    } else if (eff_num_first < d_first) {
        *finite_flag = 0.0;
        return INFINITY;
    } else {
        *finite_flag = 1.0;
        double n_val = numer[n_first];
        double d_val = denom[d_first];
        if (fabs(d_val) < 1e-15) { *finite_flag = 0.0; return NAN; }
        return n_val / d_val;
    }
}

/* ============================================================================
 * L5: Non-Unity Feedback Error Constants
 * ============================================================================ */

ErrorConstantDetail ec_nonunity_error_constants(const TransferFunction *G,
                                                  const TransferFunction *H)
{
    ErrorConstantDetail detail;
    memset(&detail, 0, sizeof(detail));
    if (!G || !H) return detail;

    double G0 = (fabs(G->denom_coeffs[0]) > 1e-15)
                ? G->gain * G->numer_coeffs[0] / G->denom_coeffs[0] : INFINITY;
    double H0 = (fabs(H->denom_coeffs[0]) > 1e-15)
                ? H->gain * H->numer_coeffs[0] / H->denom_coeffs[0] : INFINITY;

    detail.Kp_numeric = (isinf(G0) || isinf(H0)) ? INFINITY : G0 * H0;
    detail.Kp_finite = !isinf(detail.Kp_numeric);

    double lim_sG = ec_compute_Kv(G);
    double lim_sH = ec_compute_Kv(H);
    detail.Kv_numeric = lim_sG * H0 + G0 * lim_sH;
    if (fabs(detail.Kv_numeric) > 1e308) detail.Kv_numeric = INFINITY;
    detail.Kv_finite = !isinf(detail.Kv_numeric) && fabs(detail.Kv_numeric) > 1e-15;

    double lim_s2G = ec_compute_Ka(G);
    double lim_s2H = ec_compute_Ka(H);
    detail.Ka_numeric = lim_s2G * H0 + 2.0 * lim_sG * lim_sH + G0 * lim_s2H;
    if (fabs(detail.Ka_numeric) > 1e308) detail.Ka_numeric = INFINITY;
    detail.Ka_finite = !isinf(detail.Ka_numeric) && fabs(detail.Ka_numeric) > 1e-15;

    int typeG = 0, typeH = 0;
    while (typeG <= G->denom_order && fabs(G->denom_coeffs[typeG]) < 1e-15) typeG++;
    while (typeH <= H->denom_order && fabs(H->denom_coeffs[typeH]) < 1e-15) typeH++;
    detail.system_type = typeG + typeH;

    detail.Kp_inf_threshold = 1e308;
    detail.Kv_inf_threshold = 1e308;
    detail.Ka_inf_threshold = 1e308;

    return detail;
}

/* ============================================================================
 * L5: Generalized SSE, Bandwidth Estimation, Required Gain, Bode Mapping
 * ============================================================================ */

double ec_generalized_sse(const TransferFunction *E, double tol)
{
    if (!E) return NAN;

    int e_type = 0;
    while (e_type <= E->denom_order &&
           fabs(E->denom_coeffs[e_type]) < tol) e_type++;

    int n_first = 0;
    while (n_first <= E->numer_order &&
           fabs(E->numer_coeffs[n_first]) < tol) n_first++;

    int eff_num_first = n_first + 1;
    int d_first = e_type;

    if (eff_num_first > d_first) return 0.0;
    else if (eff_num_first < d_first) return INFINITY;
    else {
        double nn = E->numer_coeffs[n_first];
        double dd = E->denom_coeffs[d_first];
        if (fabs(dd) < tol) return NAN;
        return E->gain * nn / dd;
    }
}

double ec_bandwidth_error_estimate(double zeta, double omega_n,
                                    TestInputType input)
{
    if (omega_n < 1e-15) return NAN;
    if (zeta < 0.0) zeta = 0.0;

    switch (input) {
    case INPUT_STEP:
        return 0.0;
    case INPUT_RAMP:
        if (zeta < 1e-10) return 0.0;
        return 2.0 * zeta / omega_n;
    case INPUT_PARABOLA:
        return INFINITY;
    default:
        return NAN;
    }
}

double ec_required_gain(SystemType type, TestInputType input, double e_desired)
{
    if (e_desired <= 0.0) return INFINITY;

    switch (type) {
    case SYSTEM_TYPE_0:
        if (input == INPUT_STEP) return (1.0 / e_desired) - 1.0;
        return INFINITY;
    case SYSTEM_TYPE_1:
        if (input == INPUT_STEP) return 0.0;
        if (input == INPUT_RAMP) return 1.0 / e_desired;
        return INFINITY;
    case SYSTEM_TYPE_2:
        if (input == INPUT_STEP)  return 0.0;
        if (input == INPUT_RAMP)  return 0.0;
        if (input == INPUT_PARABOLA) return 1.0 / e_desired;
        return INFINITY;
    default:
        return 1.0 / e_desired;
    }
}

double ec_Kp_from_bode_low_freq(double mag_db_low)
{
    return pow(10.0, mag_db_low / 20.0);
}

double ec_Kv_from_bode_crossover(double freq_at_0db)
{
    return freq_at_0db;
}

void ec_detail_print(const ErrorConstantDetail *detail)
{
    if (!detail) return;
    printf("Error Constant Details:\n");
    printf("  System Type: %d\n", detail->system_type);
    printf("  Kp: %s", detail->Kp_finite ? "" : "INF");
    if (detail->Kp_finite) printf("%.6g", detail->Kp_numeric);
    printf("\n  Kv: %s", detail->Kv_finite ? "" : "INF");
    if (detail->Kv_finite) printf("%.6g", detail->Kv_numeric);
    printf("\n  Ka: %s", detail->Ka_finite ? "" : "INF");
    if (detail->Ka_finite) printf("%.6g", detail->Ka_numeric);
    printf("\n");
}
