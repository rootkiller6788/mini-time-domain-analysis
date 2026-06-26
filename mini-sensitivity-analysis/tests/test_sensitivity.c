/**
 * test_sensitivity.c — Comprehensive tests for mini-sensitivity-analysis
 *
 * Tests cover all core APIs across:
 *   - Complex arithmetic (L3)
 *   - Transfer function evaluation & interconnection (L2)
 *   - Sensitivity function S+T=1 identity (L1, L4)
 *   - Bode integral computation (L4)
 *   - Eigenvalue sensitivity (L4, L5)
 *   - Parameter sensitivity (L5)
 *   - Step response metrics (L6)
 *   - Second-order system sensitivity (L6)
 *   - Robustness bounds (L8)
 */

#include "sensitivity_core.h"
#include "sensitivity_transfer.h"
#include "sensitivity_parametric.h"
#include "sensitivity_bode.h"
#include "sensitivity_eigenvalue.h"
#include "sensitivity_time_domain.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define TOL 1e-6
#define TOL_SOFT 1e-2

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)
#define CHECK(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)
#define CHECK_CLOSE(a, b, tol, msg) do { if (fabs((a)-(b)) > (tol)) { FAIL(msg); return; } } while(0)

/* ==========================================================================
 * L3: Complex Number Arithmetic Tests
 * ========================================================================== */

static void test_complex_arithmetic(void) {
    TEST("complex_make");
    Complex z = complex_make(3.0, 4.0);
    CHECK(z.re == 3.0 && z.im == 4.0, "complex_make failed");
    PASS();

    TEST("complex_add");
    Complex z1 = complex_make(1.0, 2.0);
    Complex z2 = complex_make(3.0, 4.0);
    Complex sum = complex_add(z1, z2);
    CHECK_CLOSE(sum.re, 4.0, TOL, "add real");
    CHECK_CLOSE(sum.im, 6.0, TOL, "add imag");
    PASS();

    TEST("complex_sub");
    Complex diff = complex_sub(z1, z2);
    CHECK_CLOSE(diff.re, -2.0, TOL, "sub real");
    CHECK_CLOSE(diff.im, -2.0, TOL, "sub imag");
    PASS();

    TEST("complex_mul");
    Complex prod = complex_mul(z1, z2);
    CHECK_CLOSE(prod.re, -5.0, TOL, "mul real");  /* 1*3 - 2*4 */
    CHECK_CLOSE(prod.im, 10.0, TOL, "mul imag");  /* 1*4 + 2*3 */
    PASS();

    TEST("complex_div");
    Complex quot = complex_div(complex_make(5.0, 0.0), complex_make(2.0, 0.0));
    CHECK_CLOSE(quot.re, 2.5, TOL, "div real");
    CHECK_CLOSE(quot.im, 0.0, TOL, "div imag");
    /* Division by zero */
    Complex nan_div = complex_div(complex_make(1.0, 0.0), complex_make(0.0, 0.0));
    CHECK(isnan(nan_div.re) || isnan(nan_div.im), "div by zero should be NaN");
    PASS();

    TEST("complex_abs");
    CHECK_CLOSE(complex_abs(complex_make(3.0, 4.0)), 5.0, TOL, "abs 3+4j");
    CHECK_CLOSE(complex_abs(complex_make(0.0, 0.0)), 0.0, TOL, "abs 0");
    PASS();

    TEST("complex_arg");
    CHECK_CLOSE(complex_arg(complex_make(1.0, 0.0)), 0.0, TOL, "arg 1");
    CHECK_CLOSE(complex_arg(complex_make(0.0, 1.0)), M_PI/2.0, TOL, "arg j");
    CHECK_CLOSE(complex_arg(complex_make(-1.0, 0.0)), M_PI, TOL_SOFT, "arg -1");
    PASS();

    TEST("complex_conj");
    Complex conj = complex_conj(complex_make(3.0, 4.0));
    CHECK_CLOSE(conj.re, 3.0, TOL, "conj real");
    CHECK_CLOSE(conj.im, -4.0, TOL, "conj imag");
    PASS();
}

/* ==========================================================================
 * L1-L2: Sensitivity Function Tests (S+T=1 identity)
 * ========================================================================== */

static void test_sensitivity_identity(void) {
    TEST("S+T=1 identity at DC");
    /* G(s) = 1/(s+1), unity feedback */
    TransferFunction G;
    G.n = 0; G.m = 1;
    G.b[0] = 1.0;
    G.a[0] = 1.0; G.a[1] = 1.0;

    Complex L = loop_transfer(&G, NULL, complex_make(0.0, 0.0));
    Complex S = compute_sensitivity(L);
    Complex T = compute_comp_sensitivity(L);

    /* At DC: L=1, S=0.5, T=0.5, S+T=1 */
    CHECK_CLOSE(verify_st_identity(S, T), 0.0, TOL, "S+T != 1 at DC");
    PASS();

    TEST("S+T=1 identity at ω=1");
    Complex L1 = loop_transfer(&G, NULL, complex_make(0.0, 1.0));
    Complex S1 = compute_sensitivity(L1);
    Complex T1 = compute_comp_sensitivity(L1);
    CHECK_CLOSE(verify_st_identity(S1, T1), 0.0, TOL, "S+T != 1 at ω=1");
    PASS();

    TEST("S+T=1 identity at ω=100");
    Complex L100 = loop_transfer(&G, NULL, complex_make(0.0, 100.0));
    Complex S100 = compute_sensitivity(L100);
    Complex T100 = compute_comp_sensitivity(L100);
    CHECK_CLOSE(verify_st_identity(S100, T100), 0.0, TOL_SOFT, "S+T != 1 at ω=100");
    PASS();

    TEST("S for L=0 is 1");
    Complex S_zero = compute_sensitivity(complex_make(0.0, 0.0));
    CHECK_CLOSE(S_zero.re, 1.0, TOL, "S(0) != 1");
    CHECK_CLOSE(S_zero.im, 0.0, TOL, "S(0) imag != 0");
    PASS();

    TEST("T for L=0 is 0");
    Complex T_zero = compute_comp_sensitivity(complex_make(0.0, 0.0));
    CHECK_CLOSE(T_zero.re, 0.0, TOL, "T(0) != 0");
    PASS();
}

/* ==========================================================================
 * L2: Transfer Function Evaluation & Interconnection
 * ========================================================================== */

static void test_transfer_function(void) {
    TEST("tf_evaluate DC gain");
    TransferFunction G;
    G.n = 0; G.m = 1;
    G.b[0] = 5.0;
    G.a[0] = 2.0; G.a[1] = 1.0;
    Complex G0 = tf_evaluate(&G, complex_make(0.0, 0.0));
    CHECK_CLOSE(G0.re, 2.5, TOL, "DC gain wrong");  /* 5/2 */
    PASS();

    TEST("tf_evaluate at pole frequency");
    /* G(s)=1/(s+1), at s=-1 should be pole */
    TransferFunction G2;
    G2.n = 0; G2.m = 1;
    G2.b[0] = 1.0;
    G2.a[0] = 1.0; G2.a[1] = 1.0;
    Complex Gp = tf_evaluate(&G2, complex_make(-1.0, 0.0));
    /* Should return large value (pole) */
    CHECK(complex_abs(Gp) > 1e10 || isnan(Gp.re), "pole not detected");
    PASS();

    TEST("tf_series");
    TransferFunction A, B, C;
    A.n = 0; A.m = 0; A.b[0] = 2.0; A.a[0] = 1.0;
    B.n = 0; B.m = 0; B.b[0] = 3.0; B.a[0] = 1.0;
    int ret = tf_series(&A, &B, &C);
    CHECK(ret == 0, "tf_series failed");
    CHECK_CLOSE(C.b[0], 6.0, TOL, "series num wrong");  /* 2*3 = 6 */
    PASS();

    TEST("tf_parallel");
    ret = tf_parallel(&A, &B, &C);
    CHECK(ret == 0, "tf_parallel failed");
    /* G=2/1 + 3/1 = 5/1 */
    CHECK_CLOSE(C.b[0], 5.0, TOL, "parallel num wrong");
    PASS();

    TEST("tf_feedback (Black's formula)");
    /* G=A=2, H=B=3 ⇒ T = 2/(1+6) = 2/7 */
    ret = tf_feedback(&A, &B, &C);
    CHECK(ret == 0, "tf_feedback failed");
    CHECK_CLOSE(C.b[0], 2.0, TOL, "feedback num wrong");
    CHECK_CLOSE(C.a[0], 7.0, TOL, "feedback den wrong");
    PASS();

    TEST("tf_type_number (integrators)");
    TransferFunction I;
    I.n = 0; I.m = 2;
    I.b[0] = 1.0;
    I.a[0] = 0.0; I.a[1] = 0.0; I.a[2] = 1.0;
    CHECK(tf_type_number(&I) == 2, "should have 2 integrators");
    PASS();

    TEST("tf_pole_excess");
    TransferFunction PE;
    PE.n = 1; PE.m = 3;
    PE.b[0] = 1.0; PE.b[1] = 2.0;
    PE.a[0] = 1.0; PE.a[1] = 3.0; PE.a[2] = 3.0; PE.a[3] = 1.0;
    CHECK(tf_pole_excess(&PE) == 2, "pole excess should be 2");
    PASS();

    TEST("tf_dc_gain");
    CHECK_CLOSE(tf_dc_gain(&G), 2.5, TOL, "DC gain wrong");
    CHECK(isinf(tf_dc_gain(&I)), "integrator DC gain should be INF");
    PASS();
}

/* ==========================================================================
 * L2: Polynomial Arithmetic
 * ========================================================================== */

static void test_polynomial(void) {
    TEST("poly_add");
    double a[] = {1.0, 2.0, 1.0};  /* 1 + 2s + s² */
    double b[] = {3.0, 1.0};       /* 3 + s */
    double r[5];
    int deg = poly_add(a, 2, b, 1, r);
    CHECK(deg == 2, "poly_add degree wrong");
    CHECK_CLOSE(r[0], 4.0, TOL, "poly_add const");
    CHECK_CLOSE(r[1], 3.0, TOL, "poly_add s");
    CHECK_CLOSE(r[2], 1.0, TOL, "poly_add s²");
    PASS();

    TEST("poly_mul (convolution)");
    double m[5];
    deg = poly_mul(a, 2, b, 1, m);
    CHECK(deg == 3, "poly_mul degree wrong");
    /* (1+2s+s²)(3+s) = 3 + 7s + 5s² + s³ */
    CHECK_CLOSE(m[0], 3.0, TOL, "poly_mul const");
    CHECK_CLOSE(m[1], 7.0, TOL, "poly_mul s");
    CHECK_CLOSE(m[2], 5.0, TOL, "poly_mul s²");
    CHECK_CLOSE(m[3], 1.0, TOL, "poly_mul s³");
    PASS();

    TEST("poly_eval Horner");
    double p[] = {1.0, 0.0, -1.0};  /* 1 - s² */
    Complex val = poly_eval(p, 2, complex_make(2.0, 0.0));
    CHECK_CLOSE(val.re, -3.0, TOL, "poly_eval");  /* 1 - 4 = -3 */
    PASS();
}

/* ==========================================================================
 * L4-L5: Bode Integral
 * ========================================================================== */

static void test_bode_integral(void) {
    TEST("Bode integral for stable first-order system");
    /* L(s) = 2/(s+1) — stable, no RHP poles, no RHP zeros
     * Theoretical: π·ΣRe(p_k) = 0 (no RHP poles) */
    TransferFunction L;
    L.n = 0; L.m = 1;
    L.b[0] = 2.0;
    L.a[0] = 1.0; L.a[1] = 1.0;

    BodeIntegralResult result;
    bode_integral_compute(&L, 1e-4, 1e4, 201, &result);
    /* For stable system with no RHP poles, theoretical is 0 */
    CHECK_CLOSE(result.theoretical_value, 0.0, TOL, "Bode integral theoretical != 0");
    /* The numerical integral should also be near 0 */
    CHECK_CLOSE(result.integral, 0.0, 0.5, "Bode integral numerical too far from 0");
    CHECK(result.n_rhp_poles == 0, "should have 0 RHP poles");
    free(result.rhp_poles);
    free(result.rhp_zeros);
    PASS();

    TEST("Poisson kernel");
    double w = poisson_kernel(1.0, 1.0);
    CHECK_CLOSE(w, 1.0, TOL, "W(1,1) should be 2/(1+1)=1");
    double w0 = poisson_kernel(2.0, 0.0);
    CHECK_CLOSE(w0, 1.0, TOL, "W(2,0) should be 1");
    PASS();

    TEST("Sensitivity design constraints");
    SensitivityConstraints constraints;
    sensitivity_design_constraints(&L, NULL, &constraints);
    CHECK(constraints.is_feasible, "should be feasible");
    CHECK_CLOSE(constraints.min_peak_sensitivity, 1.0, TOL_SOFT, "Ms_min should be 1");
    PASS();

    TEST("Gain margin from Ms");
    double gm = gain_margin_from_Ms(2.0);
    CHECK_CLOSE(gm, 2.0, TOL_SOFT, "GM from Ms=2");  /* Ms=2 ⇒ GM≥2 */
    CHECK(isinf(gain_margin_from_Ms(1.0)), "GM should be infinite for Ms=1");
    PASS();

    TEST("Phase margin from Ms");
    double pm = phase_margin_from_Ms(2.0);
    /* PM = 2*arcsin(1/4) ≈ 0.505 rad ≈ 28.96° */
    CHECK_CLOSE(pm, 2.0 * asin(0.25), TOL, "PM from Ms=2");
    PASS();
}

/* ==========================================================================
 * L5: Parameter Sensitivity (Finite Differences)
 * ========================================================================== */

static double test_func_sq(const double *p, int n) {
    (void)n;
    return p[0] * p[0] + p[1] * p[1];  /* f = p0² + p1² */
}

static void test_fd_sensitivity(void) {
    TEST("forward FD sensitivity");
    double params[] = {3.0, 4.0};
    SensitivityResult results[2];
    fd_sensitivity_forward(test_func_sq, params, 2, NULL, results);
    CHECK_CLOSE(results[0].abs_sensitivity, 6.0, TOL_SOFT, "∂f/∂p0 wrong");  /* 2*3=6 */
    CHECK_CLOSE(results[1].abs_sensitivity, 8.0, TOL_SOFT, "∂f/∂p1 wrong");  /* 2*4=8 */
    PASS();

    TEST("central FD sensitivity (more accurate)");
    fd_sensitivity_central(test_func_sq, params, 2, NULL, results);
    CHECK_CLOSE(results[0].abs_sensitivity, 6.0, TOL, "central ∂f/∂p0 wrong");
    CHECK_CLOSE(results[1].abs_sensitivity, 8.0, TOL, "central ∂f/∂p1 wrong");
    PASS();
}

/* ==========================================================================
 * L5: Eigenvalue Computation & Sensitivity
 * ========================================================================== */

static void test_eigenvalue(void) {
    TEST("QR eigenvalues 2×2 matrix");
    double A[4] = {2.0, 1.0, 1.0, 2.0};  /* eigenvalues: 3, 1 */
    Complex evals[2];
    int ret = qr_eigenvalues(A, 2, evals);
    CHECK(ret == 0, "QR failed");

    /* Sort: check if 3 and 1 are found */
    double e1 = evals[0].re, e2 = evals[1].re;
    CHECK((fabs(e1 - 3.0) < TOL && fabs(e2 - 1.0) < TOL) ||
          (fabs(e1 - 1.0) < TOL && fabs(e2 - 3.0) < TOL),
          "2×2 eigenvalues wrong");
    PASS();

    TEST("Spectral abscissa");
    double alpha = spectral_abscissa(evals, 2);
    CHECK_CLOSE(alpha, 3.0, TOL, "spectral abscissa wrong");
    PASS();

    TEST("Spectral radius");
    double rho = spectral_radius(evals, 2);
    CHECK_CLOSE(rho, 3.0, TOL, "spectral radius wrong");
    PASS();

    TEST("Power iteration");
    double A3[9] = {3.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 1.0};
    double lambda, eigvec[3];
    int iter = power_iteration(A3, 3, 200, 1e-6, &lambda, eigvec);
    CHECK(iter > 0, "power iteration failed");
    CHECK_CLOSE(fabs(lambda), 3.0, 0.01, "dominant eigenvalue wrong");
    PASS();

    TEST("Damping ratios");
    Complex evals2[] = {
        complex_make(-1.0, 2.0),   /* ζ = 1/√5 ≈ 0.447 */
        complex_make(-1.0, -2.0)
    };
    double damp[2];
    compute_damping_ratios(evals2, 2, damp);
    CHECK_CLOSE(damp[0], 1.0/sqrt(5.0), TOL, "damping ratio wrong");
    PASS();

    TEST("Settling time estimate");
    double ts = estimate_settling_time(evals2, 2);
    CHECK_CLOSE(ts, 4.0, TOL, "settling time wrong");  /* 4/|−1| = 4 */
    PASS();
}

/* ==========================================================================
 * L6: Second-Order Step Response & Sensitivity
 * ========================================================================== */

static void test_second_order(void) {
    TEST("Second-order step response (underdamped)");
    int n = 500;
    double *t = (double *)malloc((size_t)n * sizeof(double));
    double *y = (double *)malloc((size_t)n * sizeof(double));
    CHECK(t != NULL && y != NULL, "malloc failed");

    double zeta = 0.5, omega_n = 2.0;
    second_order_step_response(zeta, omega_n, n, t, y);

    /* Check steady-state = 1 */
    CHECK_CLOSE(y[n-1], 1.0, TOL_SOFT, "steady-state not 1");

    /* Check overshoot: M_p = exp(-πζ/√(1-ζ²)) */
    double Mp_expected = exp(-M_PI * zeta / sqrt(1.0 - zeta*zeta));
    double y_max = y[0];
    for (int i = 1; i < n; i++) {
        if (y[i] > y_max) y_max = y[i];
    }
    double Mp_actual = (y_max - 1.0);  /* for unit step */
    CHECK_CLOSE(Mp_actual, Mp_expected, TOL_SOFT, "overshoot mismatch");
    PASS();

    TEST("Step response metrics");
    StepResponseMetrics metrics;
    compute_step_metrics(t, y, n, 1.0, &metrics);
    CHECK(metrics.overshoot_pct > 0.0, "overshoot should be positive");
    CHECK(metrics.rise_time > 0.0, "rise time should be positive");
    CHECK(metrics.settling_time_2pct < t[n-1], "should settle");
    PASS();

    TEST("Second-order sensitivity to ζ (analytical)");
    double *dy_dzeta = (double *)malloc((size_t)n * sizeof(double));
    CHECK(dy_dzeta != NULL, "malloc failed");
    second_order_dy_dzeta(zeta, omega_n, n, t, dy_dzeta);
    /* dy/dζ at the peak should be negative (more damping → less peak) */
    int peak_idx = 0;
    double ymax = y[0];
    for (int i = 1; i < n; i++) if (y[i] > ymax) { ymax = y[i]; peak_idx = i; }
    CHECK(dy_dzeta[peak_idx] < 0.0, "dy/dζ at peak should be negative");
    PASS();

    TEST("Second-order metric sensitivity to ω_n");
    StepMetricsSensitivity sens;
    second_order_metric_sensitivity(zeta, omega_n, &sens, 1);
    /* ∂t_p/∂ω_n = -π/(ω_n²√(1-ζ²)) < 0 */
    CHECK(sens.d_peak_time_dp < 0.0, "peak time should decrease with ω_n");
    CHECK(sens.d_settling_time_dp < 0.0, "settling time should decrease with ω_n");
    PASS();

    free(t); free(y); free(dy_dzeta);
}

/* ==========================================================================
 * L8: Robustness (Small Gain, Gap Metric)
 * ========================================================================== */

extern int small_gain_robust_stability(const TransferFunction *G_nom,
                                        const TransferFunction *K,
                                        const TransferFunction *W);
extern double nu_gap_scalar(Complex G1, Complex G2);
extern double nu_gap_metric(const TransferFunction *G1, const TransferFunction *G2,
                             double omega_min, double omega_max, int n_points);
extern double stability_margin(const TransferFunction *G, const TransferFunction *K,
                                double omega_min, double omega_max, int n_points);

static void test_robustness(void) {
    TEST("nu-gap scalar — identical plants");
    Complex g1 = complex_make(1.0, 0.0);
    double gap = nu_gap_scalar(g1, g1);
    CHECK_CLOSE(gap, 0.0, TOL, "identical plants should have 0 gap");
    PASS();

    TEST("nu-gap scalar — very different plants");
    Complex g2 = complex_make(100.0, 0.0);
    double gap2 = nu_gap_scalar(g1, g2);
    CHECK(gap2 > 0.5, "very different plants should have large gap");
    PASS();

    TEST("Small gain robust stability");
    /* G=1/(s+1), K=1, W=0.1 (sufficiently small uncertainty) */
    TransferFunction G, K, W;
    G.n = 0; G.m = 1;
    G.b[0] = 1.0;
    G.a[0] = 1.0; G.a[1] = 1.0;

    K.n = 0; K.m = 0;
    K.b[0] = 1.0; K.a[0] = 1.0;

    W.n = 0; W.m = 0;
    W.b[0] = 0.1; W.a[0] = 1.0;

    int stable = small_gain_robust_stability(&G, &K, &W);
    CHECK(stable == 1, "should be robustly stable with small uncertainty");
    PASS();
}

/* ==========================================================================
 * L1-L5: Morris Screening Method
 * ========================================================================== */

static double morris_test_func(const double *p, int n) {
    /* f = p0 + 2*p1 + 3*p2 (linear, no interaction) */
    (void)n;
    return p[0] + 2.0 * p[1] + 3.0 * p[2];
}

static void test_morris(void) {
    TEST("Morris screening — linear function");
    double bounds[6] = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
    double mu_star[3], sigma[3];
    morris_screening(morris_test_func, bounds, 3, 4, 4, mu_star, sigma);
    /* For linear function f = p0+2p1+3p2, μ* should reflect coefficients */
    CHECK(mu_star[2] > mu_star[1], "p2 should be most influential");
    CHECK(mu_star[1] > mu_star[0], "p1 should be more influential than p0");
    /* σ should be near 0 for linear function (no interactions) */
    CHECK(sigma[0] < 0.1 && sigma[1] < 0.1 && sigma[2] < 0.1,
          "sigma should be small for linear function");
    PASS();
}

/* ==========================================================================
 * L5: Normalized Sensitivity Ranking
 * ========================================================================== */

static void test_ranking(void) {
    TEST("Rank by sensitivity");
    double sens[] = {1.0, 5.0, 0.5, 3.0};
    int ranks[4];
    rank_by_sensitivity(sens, 4, ranks);
    CHECK(ranks[1] == 1, "most sensitive should be index 1");  /* 5.0 */
    CHECK(ranks[3] == 2, "second should be index 3");          /* 3.0 */
    PASS();

    TEST("Normalize sensitivities");
    double norm[4];
    normalize_sensitivities(sens, 4, norm);
    double sum = norm[0] + norm[1] + norm[2] + norm[3];
    CHECK_CLOSE(sum, 1.0, TOL, "normalized sum should be 1");
    PASS();
}

/* ==========================================================================
 * L6: State-Space Simulation
 * ========================================================================== */

static double unit_step(double t) {
    (void)t;
    return 1.0;
}

static void test_state_space(void) {
    TEST("LTI RK4 simulation");
    StateSpace ss;
    ss.n = 1; ss.m = 1; ss.p = 1;
    ss.A[0] = -1.0;  /* A[0][0] = -1.0 */
    ss.B[0] = 1.0;   /* B[0][0] = 1.0 */
    ss.C[0] = 1.0;   /* C[0][0] = 1.0 */
    ss.D[0] = 0.0;   /* D[0][0] = 0.0 */

    int n_steps = 100;
    double *t = (double *)malloc((size_t)(n_steps + 1) * sizeof(double));
    double *y = (double *)malloc((size_t)(n_steps + 1) * sizeof(double));
    double x0[] = {0.0};

    simulate_lti(&ss, 0.0, 5.0, x0, unit_step, n_steps, t, NULL, y);

    /* ẋ = -x + 1, x(0)=0 ⇒ x(t) = 1 - e^{-t}, y=x
     * At t=5, y ≈ 1 - e^{-5} ≈ 0.9933 */
    CHECK_CLOSE(y[n_steps], 1.0 - exp(-5.0), TOL_SOFT, "RK4 final value wrong");
    PASS();

    free(t); free(y);
}

/* ==========================================================================
 * L5: Mixed Sensitivity Cost
 * ========================================================================== */

static void test_mixed_sensitivity(void) {
    TEST("Mixed sensitivity cost");
    double c1 = mixed_sensitivity_cost(1.0, 0.0, 1.0, 0.0);
    CHECK_CLOSE(c1, 1.0, TOL, "|S|=1, α=1 ⇒ cost=1");
    double c2 = mixed_sensitivity_cost(0.0, 1.0, 0.0, 1.0);
    CHECK_CLOSE(c2, 1.0, TOL, "|T|=1, β=1 ⇒ cost=1");
    double c3 = mixed_sensitivity_cost(0.5, 0.5, 1.0, 1.0);
    CHECK_CLOSE(c3, 0.5, TOL, "|S|=|T|=0.5, α=β=1 ⇒ cost=0.5");
    PASS();
}

/* ==========================================================================
 * L5: Waterbed Penalty
 * ========================================================================== */

static void test_waterbed(void) {
    TEST("Waterbed penalty");
    /* If ε=0.5, ω_c=1.0, sum RHP poles=0 (stable):
     * Even stable systems have waterbed: forcing |S|≤0.5 on [0,ω_c]
     * forces |S|>1 at higher frequencies. Penalty ≈ exp(ω_c·|ln(ε)|/ω_c) = 1/ε.
     * Actually: required_above = -ω_c*ln(ε) > 0, so penalty > 1. */
    double penalty = waterbed_penalty(0.5, 1.0, 0.0);
    CHECK(penalty > 1.0, "even stable systems have waterbed penalty when |S|<1");
    /* If RHP poles present, penalty is larger */
    double penalty2 = waterbed_penalty(0.5, 1.0, 0.5);
    CHECK(penalty2 > penalty, "RHP poles increase the waterbed penalty");
    PASS();
}

/* ==========================================================================
 * Main: Run All Tests
 * ========================================================================== */

int main(void) {
    printf("=== mini-sensitivity-analysis Test Suite ===\n\n");

    printf("[L3] Complex Number Arithmetic:\n");
    test_complex_arithmetic();

    printf("\n[L1-L2] Sensitivity Function Identity:\n");
    test_sensitivity_identity();

    printf("\n[L2] Transfer Function:\n");
    test_transfer_function();

    printf("\n[L2] Polynomial Arithmetic:\n");
    test_polynomial();

    printf("\n[L4-L5] Bode Integral:\n");
    test_bode_integral();

    printf("\n[L5] Parameter Sensitivity:\n");
    test_fd_sensitivity();

    printf("\n[L5] Eigenvalue Analysis:\n");
    test_eigenvalue();

    printf("\n[L6] Second-Order Systems:\n");
    test_second_order();

    printf("\n[L8] Robustness:\n");
    test_robustness();

    printf("\n[L5] Morris Screening:\n");
    test_morris();

    printf("\n[L5] Sensitivity Ranking:\n");
    test_ranking();

    printf("\n[L6] State-Space Simulation:\n");
    test_state_space();

    printf("\n[L5] Mixed Sensitivity Cost:\n");
    test_mixed_sensitivity();

    printf("\n[L5] Waterbed Penalty:\n");
    test_waterbed();

    printf("\n========================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);

    if (tests_failed > 0) {
        printf("SOME TESTS FAILED!\n");
        return 1;
    }
    printf("All tests passed.\n");
    return 0;
}
