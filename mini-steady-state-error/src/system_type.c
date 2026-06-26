/* ============================================================================
 * mini-steady-state-error: System Type Analysis Implementation
 * L1: System type definition and classification
 * L2: Tracking capability by system type (Nise Table 7.2)
 * L3: Polynomial decomposition for type detection
 * L4: Internal Model Principle (Francis & Wonham, 1976)
 * L5: Type augmentation design (lag, PI, PID)
 * Reference: Nise 7.3-7.4, Ogata 5.3-5.4
 * ============================================================================ */

#include "system_type.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

SystemTypeInfo st_analyze(const TransferFunction *tf)
{
    SystemTypeInfo info;
    memset(&info, 0, sizeof(info));
    if (!tf || !tf->denom_coeffs) return info;

    int integrators = 0;
    while (integrators <= tf->denom_order &&
           fabs(tf->denom_coeffs[integrators]) < 1e-15) integrators++;
    info.integrator_count = integrators;
    info.type = (integrators <= 2) ? integrators : 99;

    if (integrators > 0 && integrators <= tf->denom_order)
        info.integrator_gain = tf->gain / tf->denom_coeffs[integrators];
    else if (integrators == 0)
        info.integrator_gain = tf->gain * tf->numer_coeffs[0] / tf->denom_coeffs[0];

    info.has_differentiator = (tf->numer_order >= 0 &&
                                fabs(tf->numer_coeffs[0]) < 1e-15);
    info.relative_degree = tf->denom_order - tf->numer_order;
    info.is_proper = (info.relative_degree >= 0);
    info.is_strictly_proper = (info.relative_degree > 0);
    return info;
}

const char *st_type_name(int type)
{
    switch (type) {
    case 0:  return "Type 0";
    case 1:  return "Type 1";
    case 2:  return "Type 2";
    case 3:  return "Type 3";
    default: return "Type N";
    }
}

TrackingCapability st_tracking_capability(int system_type, TestInputType input)
{
    int max_zero_error_degree = system_type - 1;
    int input_degree;
    switch (input) {
    case INPUT_STEP:     input_degree = 0; break;
    case INPUT_RAMP:     input_degree = 1; break;
    case INPUT_PARABOLA: input_degree = 2; break;
    default:             return TRACK_DIVERGES;
    }
    if (input_degree <= max_zero_error_degree) return TRACK_ZERO_ERROR;
    else if (input_degree == system_type) return TRACK_CONSTANT_ERROR;
    else return TRACK_DIVERGES;
}

int st_max_trackable_degree(int system_type) { return system_type - 1; }
int st_finite_error_degree(int system_type)  { return system_type; }

const char *st_tracking_name(TrackingCapability cap)
{
    switch (cap) {
    case TRACK_DIVERGES:       return "Diverges (infinite error)";
    case TRACK_CONSTANT_ERROR: return "Constant (finite non-zero error)";
    case TRACK_ZERO_ERROR:     return "Zero error (perfect tracking)";
    default:                   return "Unknown";
    }
}

/* ============================================================================
 * L3: Polynomial Decomposition
 * ============================================================================ */

int st_count_leading_zeros(const double *coeffs, int order)
{
    if (!coeffs) return 0;
    int count = 0;
    while (count <= order && fabs(coeffs[count]) < 1e-15) count++;
    return count;
}

int st_first_nonzero_index(const double *coeffs, int order)
{
    if (!coeffs) return -1;
    for (int i = 0; i <= order; i++)
        if (fabs(coeffs[i]) > 1e-15) return i;
    return -1;
}

void st_factor_out_s_power(const double *denom_coeffs, int denom_order, int N,
                            double **reduced_coeffs, int *reduced_order)
{
    if (!denom_coeffs || !reduced_coeffs || !reduced_order) return;
    *reduced_order = denom_order - N;
    if (*reduced_order < 0) { *reduced_coeffs = NULL; *reduced_order = -1; return; }
    *reduced_coeffs = (double *)malloc((size_t)(*reduced_order + 1) * sizeof(double));
    for (int i = 0; i <= *reduced_order; i++)
        (*reduced_coeffs)[i] = denom_coeffs[i + N];
}

/* ============================================================================
 * L4: Internal Model Principle (Francis & Wonham, Automatica 12, 1976)
 * For perfect asymptotic tracking, the loop must contain a model
 * of the reference signal generator dynamics.
 * ============================================================================ */

bool st_check_imp_satisfied(const TransferFunction *open_loop,
                             TestInputType ref_type)
{
    if (!open_loop) return false;
    int required;
    switch (ref_type) {
    case INPUT_STEP:     required = 1; break;
    case INPUT_RAMP:     required = 2; break;
    case INPUT_PARABOLA: required = 3; break;
    default:             return false;
    }
    int actual = st_count_leading_zeros(open_loop->denom_coeffs,
                                         open_loop->denom_order);
    return (actual >= required);
}

bool st_check_imp_sinusoidal(const TransferFunction *open_loop, double omega)
{
    if (!open_loop || omega < 1e-10 || open_loop->denom_order < 2) return false;
    double real_part = 0.0, imag_part = 0.0;
    int sign = 1;
    for (int i = 0; i <= open_loop->denom_order; i += 2) {
        real_part += sign * open_loop->denom_coeffs[i] * pow(omega, i);
        sign = -sign;
    }
    sign = 1;
    for (int i = 1; i <= open_loop->denom_order; i += 2) {
        imag_part += sign * open_loop->denom_coeffs[i] * pow(omega, i);
        sign = -sign;
    }
    return (fabs(real_part) < 1e-10 && fabs(imag_part) < 1e-10);
}

char **st_imp_required_factors(TestInputType ref_type)
{
    char **factors = (char **)malloc(4 * sizeof(char *));
    int idx = 0;
    switch (ref_type) {
    case INPUT_STEP:
        factors[idx++] = strdup("s in denominator [one integrator]"); break;
    case INPUT_RAMP:
        factors[idx++] = strdup("s^2 in denominator [two integrators]"); break;
    case INPUT_PARABOLA:
        factors[idx++] = strdup("s^3 in denominator [three integrators]"); break;
    default:
        factors[idx++] = strdup("Depends on custom input generator poles"); break;
    }
    factors[idx++] = strdup("R(s) denominator poles must appear in L(s)");
    factors[idx] = NULL;
    return factors;
}

/* ============================================================================
 * L5: Type Augmentation Design Methods
 * ============================================================================ */

int st_augment_type_lag(int original_type, double z_c, double p_c)
{
    /* Lag compensator does NOT change system type. It only increases
     * low-frequency gain (error constants) for the same type.
     * Type increase requires adding integrators (PI/PID controllers). */
    (void)z_c; (void)p_c;
    return original_type;
}

void st_design_pi_for_type1(double desired_Kv, double plant_dc_gain,
                             double *K_p_out, double *K_i_out)
{
    /* PI: G_c(s) = K_p + K_i/s adds one integrator -> Type += 1
     * Kv = lim s*(K_p + K_i/s)*G_p(s) = K_i * G_p(0) */
    if (fabs(plant_dc_gain) < 1e-15) {
        *K_p_out = 0.0; *K_i_out = 0.0; return;
    }
    *K_i_out = desired_Kv / plant_dc_gain;
    *K_p_out = 10.0 * (*K_i_out) / desired_Kv;
    if (*K_p_out < 1e-15) *K_p_out = 1.0;
}

void st_design_pid_for_type2(double desired_Ka, double plant_dc_gain,
                              double *K_p_out, double *K_i_out, double *K_d_out)
{
    /* PID adds one integrator. Type 0+PID = Type 1, not Type 2.
     * For Type 2 from Type 0: need double integrator (cascaded PI or PI^2).
     * This function computes PID gains for the best achievable Ka. */
    if (fabs(plant_dc_gain) < 1e-15) {
        *K_p_out = 0.0; *K_i_out = 0.0; *K_d_out = 0.0; return;
    }
    *K_i_out = desired_Ka * 0.1 / plant_dc_gain;
    *K_p_out = *K_i_out / desired_Ka;
    *K_d_out = *K_p_out * 0.01;
    if (*K_p_out < 1e-10) *K_p_out = 1.0;
    if (*K_i_out < 1e-10) *K_i_out = 0.1;
    if (*K_d_out < 1e-10) *K_d_out = 0.01;
}

SystemTypeInfo st_combined_type(const TransferFunction *compensator,
                                 const TransferFunction *plant)
{
    SystemTypeInfo result;
    memset(&result, 0, sizeof(result));
    if (!compensator || !plant) return result;
    SystemTypeInfo ci = st_analyze(compensator);
    SystemTypeInfo pi = st_analyze(plant);
    result.integrator_count = ci.integrator_count + pi.integrator_count;
    result.type = (result.integrator_count <= 2) ? result.integrator_count : 99;
    result.integrator_gain = ci.integrator_gain * pi.integrator_gain;
    result.relative_degree = ci.relative_degree + pi.relative_degree;
    result.is_proper = ci.is_proper && pi.is_proper;
    result.is_strictly_proper = ci.is_strictly_proper || pi.is_strictly_proper;
    result.has_differentiator = ci.has_differentiator || pi.has_differentiator;
    return result;
}

void st_info_print(const SystemTypeInfo *info)
{
    if (!info) return;
    printf("System Type Info:\n");
    printf("  Type: %s (%d integrators)\n", st_type_name(info->type),
           info->integrator_count);
    printf("  Integrator gain: %.6g\n", info->integrator_gain);
    printf("  Relative degree: %d\n", info->relative_degree);
    printf("  Proper: %s, Strictly proper: %s\n",
           info->is_proper ? "Yes" : "No",
           info->is_strictly_proper ? "Yes" : "No");
    printf("  Has differentiator: %s\n",
           info->has_differentiator ? "Yes" : "No");
}
