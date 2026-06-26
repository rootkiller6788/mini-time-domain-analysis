/* ============================================================================
 * mini-steady-state-error: Sensitivity Analysis Implementation
 *
 * L2: Sensitivity function S(s) = 1/(1+L(s)), complementary sensitivity T(s)
 * L3: Frequency-domain sensitivity evaluation
 * L5: Parameter sensitivity, condition number analysis
 * L6: Bode integral (waterbed effect), Bode-Freudenberg-Looze theorem
 * L8: Structured singular value (mu) for robust SSE analysis
 *
 * Reference: Nise 7.7, Ogata 5.8, Skogestad & Postlethwaite (2005)
 * ============================================================================ */

#include "sensitivity_analysis.h"
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
 * L2: Sensitivity Function Computation
 *
 * S(s) = 1/(1 + L(s))  where L(s) = G(s)H(s)
 * T(s) = L(s)/(1 + L(s)) = 1 - S(s) (for unity feedback)
 *
 * For polynomial form: S(s) = D_L(s) / (D_L(s) + N_L(s))
 *                     T(s) = N_L(s) / (D_L(s) + N_L(s))
 * ============================================================================ */

TransferFunction sa_compute_sensitivity(const TransferFunction *open_loop)
{
    TransferFunction result = {NULL, -1, NULL, -1, 0.0};
    if (!open_loop) return result;

    /* S(s) = 1/(1+L(s)) = D_L/(D_L + N_L)
     * Denominator order = max(denom_order, numer_order) */
    int max_order = (open_loop->numer_order > open_loop->denom_order)
                     ? open_loop->numer_order : open_loop->denom_order;
    result = tf_alloc(max_order, max_order);
    if (result.numer_order < 0) return result;

    result.gain = 1.0;
    int i;

    /* Numerator = D_L(s) */
    for (i = 0; i <= open_loop->denom_order; i++)
        result.numer_coeffs[i] = open_loop->denom_coeffs[i];

    /* Denominator = D_L(s) + L_gain * N_L(s) */
    for (i = 0; i <= open_loop->denom_order; i++)
        result.denom_coeffs[i] = open_loop->denom_coeffs[i];
    for (i = 0; i <= open_loop->numer_order; i++)
        result.denom_coeffs[i] += open_loop->gain * open_loop->numer_coeffs[i];

    return result;
}

TransferFunction sa_compute_complementary_sensitivity(const TransferFunction *open_loop)
{
    TransferFunction result = {NULL, -1, NULL, -1, 0.0};
    if (!open_loop) return result;

    /* T(s) = L/(1+L) = N_L/(D_L + N_L) */
    int max_order = (open_loop->numer_order > open_loop->denom_order)
                     ? open_loop->numer_order : open_loop->denom_order;
    result = tf_alloc(open_loop->numer_order, max_order);
    if (result.numer_order < 0) return result;

    result.gain = 1.0;
    int i;

    /* Numerator = open_loop->gain * N_L(s) */
    for (i = 0; i <= open_loop->numer_order; i++)
        result.numer_coeffs[i] = open_loop->gain * open_loop->numer_coeffs[i];

    /* Denominator = D_L(s) + open_loop->gain * N_L(s) */
    for (i = 0; i <= open_loop->denom_order; i++)
        result.denom_coeffs[i] = open_loop->denom_coeffs[i];
    for (i = 0; i <= open_loop->numer_order; i++)
        result.denom_coeffs[i] += open_loop->gain * open_loop->numer_coeffs[i];

    return result;
}

double sa_dc_sensitivity(const ErrorConstants *ec)
{
    if (!ec) return NAN;
    /* S(0) = 1/(1+Kp) for Type 0, S(0) = 0 for Type >= 1 */
    if (ec->Kp_is_inf) return 0.0;
    return 1.0 / (1.0 + ec->Kp);
}

/* ============================================================================
 * L3: Frequency-Domain Sensitivity at omega
 *
 * Evaluate |S(jw)| = 1 / |1 + L(jw)|
 * L(jw) = K * N(jw) / D(jw)
 * ============================================================================ */

double sa_sensitivity_at_freq(const TransferFunction *open_loop, double omega)
{
    if (!open_loop) return NAN;

    /* Evaluate N(jw) and D(jw):
     * N(jw) = n0 + n1*(jw) + n2*(jw)^2 + ...
     * Separate real and imaginary parts: even powers are real, odd are imag */
    double Nr = 0.0, Ni = 0.0, Dr = 0.0, Di = 0.0;
    int i, sign = 1;

    for (i = 0; i <= open_loop->numer_order; i++) {
        if (i % 2 == 0) Nr += sign * open_loop->numer_coeffs[i] * pow(omega, i);
        else            Ni += sign * open_loop->numer_coeffs[i] * pow(omega, i);
        sign = -sign; /* j^(i) cycling: even->real, odd->imag with alternating signs */
    }

    sign = 1;
    for (i = 0; i <= open_loop->denom_order; i++) {
        if (i % 2 == 0) Dr += sign * open_loop->denom_coeffs[i] * pow(omega, i);
        else            Di += sign * open_loop->denom_coeffs[i] * pow(omega, i);
        sign = -sign;
    }

    /* L(jw) = K * (Nr + j*Ni) / (Dr + j*Di) */
    double Lr = open_loop->gain * (Nr * Dr + Ni * Di) / (Dr * Dr + Di * Di);
    double Li = open_loop->gain * (Ni * Dr - Nr * Di) / (Dr * Dr + Di * Di);

    /* S(jw) = 1/(1+L(jw)) = 1/((1+Lr) + j*Li)
     * |S| = 1 / sqrt((1+Lr)^2 + Li^2) */
    double denom = sqrt((1.0 + Lr) * (1.0 + Lr) + Li * Li);
    if (denom < 1e-15) return INFINITY;
    return 1.0 / denom;
}

/* ============================================================================
 * L5: Parameter Sensitivity Analysis
 *
 * Sensitivity of steady-state error to parameter variations.
 * For a parameter p and error e_ss:
 *   Sensitivity S_p = (de_ss/dp) * (p/e_ss)  [normalized]
 * or simply de_ss/dp [absolute].
 *
 * For Type 0 + step: e_ss = 1/(1+K), de_ss/dK = -1/(1+K)^2 = -e_ss^2
 * For Type 1 + ramp: e_ss = 1/K_v, de_ss/dK_v = -1/K_v^2 = -e_ss/K_v
 * ============================================================================ */

double sa_gain_sensitivity(int system_type, TestInputType input_type,
                            double K, double e_ss_current)
{
    if (e_ss_current <= 0.0 || isinf(e_ss_current)) return 0.0;

    switch (system_type) {
    case 0:
        if (input_type == INPUT_STEP)
            return -e_ss_current * e_ss_current;
        return 0.0;
    case 1:
        if (input_type == INPUT_RAMP && fabs(K) > 1e-15)
            return -e_ss_current / K;
        return 0.0;
    case 2:
        if (input_type == INPUT_PARABOLA && fabs(K) > 1e-15)
            return -e_ss_current / K;
        return 0.0;
    default:
        return 0.0;
    }
}

ParameterSensitivity sa_parameter_sensitivity(const TransferFunction *G_nominal,
                                                double param_delta)
{
    ParameterSensitivity ps;
    memset(&ps, 0, sizeof(ps));

    if (!G_nominal || param_delta <= 0.0) return ps;

    /* Nominal SSE */
    ErrorConstants ec_nom = sse_error_constants_from_tf(G_nominal);
    double e_nom = ec_nom.Kp_is_inf ? 0.0 : 1.0/(1.0 + ec_nom.Kp);

    /* Perturb gain by param_delta */
    TransferFunction G_pert = tf_copy(G_nominal);
    if (G_pert.numer_order < 0) return ps;

    G_pert.gain *= (1.0 + param_delta);
    ErrorConstants ec_pert = sse_error_constants_from_tf(&G_pert);
    double e_pert = ec_pert.Kp_is_inf ? 0.0 : 1.0/(1.0 + ec_pert.Kp);

    if (fabs(param_delta) > 1e-15)
        ps.d_ess_d_K = (e_pert - e_nom) / (G_nominal->gain * param_delta);
    else
        ps.d_ess_d_K = 0.0;

    /* Perturb dominant pole (denom_coeffs[1]) */
    if (G_nominal->denom_order >= 1 && fabs(G_nominal->denom_coeffs[1]) > 1e-15) {
        double orig_d1 = G_nominal->denom_coeffs[1];
        G_pert.denom_coeffs[1] *= (1.0 + param_delta);
        ec_pert = sse_error_constants_from_tf(&G_pert);
        e_pert = ec_pert.Kp_is_inf ? 0.0 : 1.0/(1.0 + ec_pert.Kp);
        ps.d_ess_d_pole = (e_pert - e_nom) / (orig_d1 * param_delta);
        G_pert.denom_coeffs[1] = orig_d1; /* restore */
    }

    /* Perturb zero location (numer_coeffs[1] if exists) */
    if (G_nominal->numer_order >= 1 && fabs(G_nominal->numer_coeffs[1]) > 1e-15) {
        double orig_n1 = G_nominal->numer_coeffs[1];
        G_pert.numer_coeffs[1] *= (1.0 + param_delta);
        ec_pert = sse_error_constants_from_tf(&G_pert);
        e_pert = ec_pert.Kp_is_inf ? 0.0 : 1.0/(1.0 + ec_pert.Kp);
        ps.d_ess_d_zero = (e_pert - e_nom) / (orig_n1 * param_delta);
    }

    /* Worst-case within perturbation range */
    ps.worst_case_e_ss = e_nom + fabs(param_delta) *
        (fabs(ps.d_ess_d_K) + fabs(ps.d_ess_d_pole) + fabs(ps.d_ess_d_zero));

    /* Variance-based normalized sensitivity */
    ps.param_variance = (fabs(ps.d_ess_d_K) + fabs(ps.d_ess_d_pole)
                          + fabs(ps.d_ess_d_zero)) * param_delta / (e_nom + 1e-10);

    tf_free(&G_pert);
    return ps;
}

double sa_error_condition_number(double e_ss, double d_ess_d_param,
                                  double param_value)
{
    if (fabs(e_ss) < 1e-15) return INFINITY;
    return fabs(d_ess_d_param * param_value / e_ss);
}

/* ============================================================================
 * L6: Bode Sensitivity Integral (Waterbed Effect)
 *
 * Bode (1945): For stable open-loop L(s) with relative degree >= 2:
 *   integral_0^inf ln|S(jw)| dw = 0
 *
 * For unstable plants with RHP poles p_i:
 *   integral_0^inf ln|S(jw)| dw = pi * sum_i Re(p_i)
 *
 * Freudenberg & Looze (1985): For non-minimum-phase systems with
 *   RHP zeros z_j, there are additional constraints.
 *   integral_0^inf ln|S(jw)| * W(z_j, w) dw = pi * ln|B_p^{-1}(z_j)|
 *
 * This function approximates the integral via trapezoidal quadrature
 * over a finite frequency range [omega_min, omega_max].
 * ============================================================================ */

BodeIntegral sa_bode_integral(const TransferFunction *open_loop,
                               const double *rhp_poles, int num_rhp,
                               double omega_min, double omega_max,
                               int num_points)
{
    BodeIntegral bi;
    memset(&bi, 0, sizeof(bi));

    if (!open_loop || num_points < 2) return bi;

    /* Expected integral from RHP poles */
    double expected = 0.0;
    int i;
    for (i = 0; i < num_rhp; i++) {
        if (rhp_poles[i] > 0) expected += rhp_poles[i];
    }
    expected *= M_PI;

    /* Numerical integration of ln|S(jw)| */
    double dw = (omega_max - omega_min) / (num_points - 1);
    double integral = 0.0;
    double pos_area = 0.0, neg_area = 0.0;
    double prev_lnS = 0.0;
    bool first = true;

    for (i = 0; i < num_points; i++) {
        double w = omega_min + i * dw;
        if (w < 1e-15) w = 1e-10; /* avoid singularity at w=0 */
        double S_mag = sa_sensitivity_at_freq(open_loop, w);
        if (S_mag < 1e-15) S_mag = 1e-10;
        if (isinf(S_mag)) S_mag = 1e10;

        double lnS = log(S_mag);
        if (lnS > 0) pos_area += lnS * dw;
        else         neg_area += lnS * dw;

        if (!first) {
            integral += 0.5 * (prev_lnS + lnS) * dw;
        }
        prev_lnS = lnS;
        first = false;
    }

    bi.waterbed_integral = integral;
    bi.positive_area = pos_area;
    bi.negative_area = neg_area;
    bi.waterbed_violated = (num_rhp > 0) && (fabs(integral - expected) > 0.1 * expected);

    return bi;
}

double sa_waterbed_peak_bound(double reduction_db, double bandwidth,
                               int relative_degree)
{
    (void)bandwidth;
    /* If |S| < M_low (< 1) from 0 to w_a, then |S| must exceed M_high > 1
     * at some frequency above w_a.
     *
     * Bode integral gives constraint: w_a * |ln(M_low)| <= (w_b - w_a) * ln(M_high)
     * where w_b is some frequency above w_a.
     *
     * Simplification (Skogestad & Postlethwaite, p.62):
     * M_high >= exp(w_a * |ln(M_low)| / (w_b - w_a))
     *
     * Using w_b = w_a * 2 (typical): M_high >= exp(|ln(M_low)|) = 1/M_low
     *
     * More generally: M_peak >= exp(pi * relative_degree * |reduction_dB| / 40) */

    double reduction_linear = pow(10.0, -fabs(reduction_db) / 20.0);
    double M_low = (reduction_linear < 1.0) ? reduction_linear : 1.0/reduction_linear;

    /* Waterbed bound with relative degree consideration */
    double bound = exp(fabs(log(M_low)) * (1.0 + 0.1 * relative_degree));

    /* Convert back to dB */
    return 20.0 * log10(bound);
}

/* ============================================================================
 * L8: Structured Singular Value (mu) Analysis
 *
 * The structured singular value mu(M) is defined as:
 *   mu(M) = 1 / min{sigma_bar(Delta) : det(I - M*Delta) = 0, Delta in structure}
 *
 * For robust steady-state error, mu quantifies the smallest structured
 * uncertainty that can violate error specifications.
 *
 * This implementation computes an upper bound for mu using the Osborne
 * iteration (scaled singular value).
 * ============================================================================ */

static double matrix_max_singular_value(const double *A, int n)
{
    /* Power iteration for maximum singular value of A.
     * sigma_max = sqrt(max eigenvalue of A^T * A) */
    double *v = (double *)calloc((size_t)n, sizeof(double));
    double *Av = (double *)calloc((size_t)n, sizeof(double));
    double *ATv = (double *)calloc((size_t)n, sizeof(double));
    int i, j;

    /* Initialize random vector */
    for (i = 0; i < n; i++) v[i] = 1.0 / sqrt((double)n);

    double sigma = 0.0, prev_sigma = 0.0;
    int iter;
    for (iter = 0; iter < 100; iter++) {
        /* Av = A * v */
        for (i = 0; i < n; i++) {
            Av[i] = 0.0;
            for (j = 0; j < n; j++) Av[i] += A[i * n + j] * v[j];
        }
        /* ATv = A^T * Av */
        for (i = 0; i < n; i++) {
            ATv[i] = 0.0;
            for (j = 0; j < n; j++) ATv[i] += A[j * n + i] * Av[j];
        }
        /* New sigma estimate */
        double norm = 0.0;
        for (i = 0; i < n; i++) norm += ATv[i] * ATv[i];
        norm = sqrt(norm);
        if (norm < 1e-15) break;
        sigma = sqrt(norm); /* sqrt of Rayleigh quotient */

        if (fabs(sigma - prev_sigma) < 1e-10) break;
        prev_sigma = sigma;

        /* v = ATv / norm */
        for (i = 0; i < n; i++) v[i] = ATv[i] / norm;
    }

    free(v); free(Av); free(ATv);
    return sigma;
}

static double matrix_inf_norm(const double *A, int n)
{
    double max_row = 0.0;
    int i, j;
    for (i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (j = 0; j < n; j++) row_sum += fabs(A[i * n + j]);
        if (row_sum > max_row) max_row = row_sum;
    }
    return max_row;
}

StructuredSSV sa_structured_sv(const double *M, int dim,
                                const int *block_structure, int num_blocks,
                                int max_iter, double tol)
{
    StructuredSSV result;
    memset(&result, 0, sizeof(result));

    if (!M || dim <= 0 || !block_structure || num_blocks <= 0) return result;

    /* Upper bound via scaled singular value (Osborne iteration):
     * mu_upper = inf_{D in D_structure} sigma_max(D*M*D^{-1})
     *
     * D = diag(d1*I_{k1}, d2*I_{k2}, ..., dN*I_{kN})  -- block-diagonal scaling
     *
     * Simplified: use Perron-Frobenius eigenvalue upper bound.
     * For full-block uncertainty: mu(M) <= sigma_max(M).
     * For repeated scalar: mu(M) <= spectral_radius(M). */

    double *M_copy = (double *)malloc((size_t)(dim * dim) * sizeof(double));
    memcpy(M_copy, M, (size_t)(dim * dim) * sizeof(double));

    /* Start with identity scaling */
    double *D = (double *)calloc((size_t)dim, sizeof(double));
    int i;
    for (i = 0; i < dim; i++) D[i] = 1.0;

    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        /* Osborne iteration: balance the matrix */
        double max_change = 0.0;
        int blk_start = 0;
        int b;
        for (b = 0; b < num_blocks; b++) {
            int blk_size = block_structure[b];
            /* Sum of absolute values in block row and column */
            double row_sum = 0.0, col_sum = 0.0;
            int r, c;
            for (r = blk_start; r < blk_start + blk_size; r++) {
                for (c = 0; c < dim; c++) {
                    if (c < blk_start || c >= blk_start + blk_size) {
                        row_sum += fabs(M_copy[r * dim + c]);
                        col_sum += fabs(M_copy[c * dim + r]);
                    }
                }
            }
            if (row_sum > 1e-15 && col_sum > 1e-15) {
                double f = sqrt(row_sum / col_sum);
                if (fabs(f - 1.0) > max_change) max_change = fabs(f - 1.0);
                /* Scale block */
                for (r = blk_start; r < blk_start + blk_size; r++) {
                    int j;
                    for (j = 0; j < dim; j++) {
                        M_copy[r * dim + j] *= f;
                        M_copy[j * dim + r] /= f;
                    }
                }
                D[blk_start] *= f;
            }
            blk_start += blk_size;
        }
        if (max_change < tol) break;
    }

    double mu_ub = matrix_max_singular_value(M_copy, dim);
    double mu_lb = matrix_inf_norm(M_copy, dim) / (double)dim;

    result.mu_upper_bound = mu_ub;
    result.mu_lower_bound = (mu_lb < mu_ub) ? mu_lb : 0.0;
    result.robust_margin = (mu_ub > 1e-15) ? 1.0 / mu_ub : INFINITY;
    result.perturbation_count = num_blocks;

    /* Nominal SSE estimate (DC gain of T(s)) */
    result.nominal_e_ss = 0.0; /* depends on system type */
    result.worst_case_e_ss_mu = result.nominal_e_ss + 1.0 / result.robust_margin;

    free(M_copy); free(D);
    return result;
}

double sa_robust_sse_bound(const TransferFunction *G_nominal,
                            const TransferFunction *W_uncertainty,
                            const double omega_range[2])
{
    if (!G_nominal || !W_uncertainty || !omega_range) return NAN;

    /* Robust SSE bound: max SSE over all allowed perturbations.
     * Using the Small Gain Theorem: if ||W*Delta||_inf < 1,
     * the worst-case sensitivity peak is bounded.
     *
     * At DC: worst-case e_ss = |S_nom(0)| / (1 - |W(0)*T_nom(0)|)
     *
     * Simplified: e_ss_wc = e_ss_nom / (1 - |W(0)|) */

    double S0 = 1.0 / (1.0 + sse_dc_gain(G_nominal));
    double W0 = sse_dc_gain(W_uncertainty);

    if (fabs(W0) >= 1.0 - 1e-10) return INFINITY; /* robustness lost */
    return S0 / (1.0 - fabs(W0));
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

void sa_metrics_print(const SensitivityMetrics *m)
{
    if (!m) return;
    printf("Sensitivity Metrics:\n");
    printf("  S_dc: %.6g (%.2f dB)\n", m->S_dc, 20.0*log10(m->S_dc));
    printf("  S_peak (Ms): %.6g\n", m->S_peak);
    printf("  T_dc: %.6g\n", m->T_dc);
    printf("  T_peak (Mt): %.6g\n", m->T_peak);
    printf("  Bandwidth: %.4g rad/s\n", m->S_bandwidth);
}

void sa_paramsens_print(const ParameterSensitivity *ps)
{
    if (!ps) return;
    printf("Parameter Sensitivity:\n");
    printf("  dE_ss/dK:   %.6g\n", ps->d_ess_d_K);
    printf("  dE_ss/dPole: %.6g\n", ps->d_ess_d_pole);
    printf("  dE_ss/dZero: %.6g\n", ps->d_ess_d_zero);
    printf("  Worst-case e_ss: %.6g\n", ps->worst_case_e_ss);
    printf("  Variance sensitivity: %.6g\n", ps->param_variance);
}

void sa_bode_print(const BodeIntegral *bi)
{
    if (!bi) return;
    printf("Bode Integral:\n");
    printf("  Integral: %.6g\n", bi->waterbed_integral);
    printf("  Positive area (|S|>1): %.6g\n", bi->positive_area);
    printf("  Negative area (|S|<1): %.6g\n", bi->negative_area);
    printf("  Waterbed violated: %s\n", bi->waterbed_violated ? "YES" : "NO");
}
