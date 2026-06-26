/**
 * @file test_impulse_step.c
 * @brief Comprehensive tests for impulse/step response module.
 *
 * Tests cover:
 *   L1: Core definitions (struct sizes, damping classification)
 *   L2: System model utilities (pole computation, stability)
 *   L4: Mathematical theorems (impulse-step relationship, Final Value Theorem)
 *   L5: Response computation accuracy (closed-form vs numeric)
 *   L6: Canonical problems (first-order, second-order step/impulse)
 *   L7: Applications (DC motor, vehicle, thermal)
 */

#include "time_response_common.h"
#include "impulse_response.h"
#include "step_response.h"
#include "convolution_response.h"
#include "response_metrics.h"
#include "discrete_response.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST: %-50s ... ", name); } while(0)
#define PASS()     do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg)  do { printf("FAIL: %s\n", msg); } while(0)

#define TOL 1e-6
#define TOL_L 1e-3  /* Loose tolerance for numerical integration */

/* ==========================================================================
 * L1 Tests: Core Definitions and Damping Classification
 * ========================================================================== */

static void test_damping_classification(void) {
    TEST("Damping classification");

    SecondOrderModel undamped      = {1.0, 0.0, 1.0};
    SecondOrderModel underdamped   = {1.0, 0.3, 1.0};
    SecondOrderModel critical      = {1.0, 1.0, 1.0};
    SecondOrderModel overdamped    = {1.0, 2.0, 1.0};
    SecondOrderModel unstable      = {1.0, -0.1, 1.0};

    assert(classify_damping(&undamped)    == DAMPING_UNDAMPED);
    assert(classify_damping(&underdamped) == DAMPING_UNDERDAMPED);
    assert(classify_damping(&critical)    == DAMPING_CRITICALLY_DAMPED);
    assert(classify_damping(&overdamped)  == DAMPING_OVERDAMPED);
    assert(classify_damping(&unstable)    == DAMPING_UNSTABLE);

    PASS();
}

static void test_second_order_poles(void) {
    TEST("Second-order pole computation");

    SecondOrderModel sys = {1.0, 0.5, 2.0};  /* zeta=0.5, wn=2 */
    double s1[2], s2[2];
    second_order_poles(&sys, s1, s2);

    /* Expected: s = -1 +/- j*sqrt(3) = -1 +/- 1.73205j */
    assert(fabs(s1[0] - (-1.0)) < TOL);
    assert(fabs(s1[1] - (sqrt(3.0))) < TOL);
    assert(fabs(s2[0] - (-1.0)) < TOL);
    assert(fabs(s2[1] - (-sqrt(3.0))) < TOL);

    /* Overdamped: zeta=2, wn=1 => s = -2 +/- sqrt(3) */
    SecondOrderModel ov = {1.0, 2.0, 1.0};
    second_order_poles(&ov, s1, s2);
    assert(fabs(s2[0] - (-2.0 - sqrt(3.0))) < TOL);  /* more negative pole */
    assert(fabs(s2[1]) < TOL);
    assert(fabs(s1[0] - (-2.0 + sqrt(3.0))) < TOL);  /* less negative pole */
    assert(fabs(s1[1]) < TOL);

    PASS();
}

static void test_design_formulas(void) {
    TEST("Design formulas (damping <-> overshoot)");

    /* 20% overshoot -> zeta ≈ 0.456 (Ogata Table 5-2) */
    double zeta = damping_from_overshoot(20.0);
    assert(fabs(zeta - 0.456) < 0.01);

    /* No overshoot -> zeta >= 1.0 */
    assert(damping_from_overshoot(0.0) >= 1.0);

    /* omega_n from settling time: ts=2s, zeta=0.5 -> omega_n=4/(0.5*2)=4 */
    double wn = omega_n_from_settling_time(2.0, 0.5);
    assert(fabs(wn - 4.0) < TOL);

    /* DC gain from steady state */
    assert(fabs(dc_gain_from_steady_state(5.0, 1.0) - 5.0) < TOL);

    PASS();
}

/* ==========================================================================
 * L4 Tests: Mathematical Theorems
 * ========================================================================== */

static void test_impulse_step_relationship(void) {
    TEST("Impulse-step integral/differential relationship");

    /* For first-order: d/dt(step) = impulse */
    FirstOrderModel fo = {2.0, 0.5};
    double t = 0.3;

    /* Numerical derivative of step at t */
    double dt_small = 1e-6;
    double step_plus = first_order_step(&fo, t + dt_small);
    double step_minus = first_order_step(&fo, t - dt_small);
    double num_deriv = (step_plus - step_minus) / (2.0 * dt_small);

    double impulse_val = first_order_impulse(&fo, t);
    assert(fabs(num_deriv - impulse_val) < 1e-3);

    PASS();
}

static void test_final_value_theorem(void) {
    TEST("Final Value Theorem");

    /* FVT: y_ss = lim_{t->inf} y_step(t) = lim_{s->0} G(s) = K */
    SecondOrderModel so = {3.0, 0.7, 2.0};
    assert(fabs(second_order_steady_state(&so) - 3.0) < TOL);

    /* Verify with large t */
    double y_large_t = second_order_step(&so, 100.0);
    assert(fabs(y_large_t - 3.0) < 1e-4);

    PASS();
}

static void test_first_order_63_percent(void) {
    TEST("First-order 63.2% at t=tau");

    FirstOrderModel fo = {1.0, 1.0};
    double y_at_tau = first_order_step(&fo, 1.0);
    double expected = 1.0 - 1.0 / M_E;
    assert(fabs(y_at_tau - expected) < TOL);

    PASS();
}

static void test_stability_check(void) {
    TEST("Routh-Hurwitz stability check");

    /* Stable first-order: G(s) = 1/(s+1), den=[1,1] */
    TransferFunction tf_stable;
    double num1[1] = {1.0};
    double den1[2] = {1.0, 1.0};
    tf_stable.num = num1; tf_stable.num_deg = 0;
    tf_stable.den = den1; tf_stable.den_deg = 1;
    assert(transfer_function_is_stable(&tf_stable) == 1);

    /* Unstable: G(s) = 1/(s-1), den=[-1,1] */
    TransferFunction tf_unstable;
    double den2[2] = {-1.0, 1.0};
    tf_unstable.num = num1; tf_unstable.num_deg = 0;
    tf_unstable.den = den2; tf_unstable.den_deg = 1;
    assert(transfer_function_is_stable(&tf_unstable) == 0);

    /* Stable second-order: s^2+2s+1, den=[1,2,1] */
    TransferFunction tf_so;
    double num3[1] = {1.0};
    double den3[3] = {1.0, 2.0, 1.0};
    tf_so.num = num3; tf_so.num_deg = 0;
    tf_so.den = den3; tf_so.den_deg = 2;
    assert(transfer_function_is_stable(&tf_so) == 1);

    PASS();
}

/* ==========================================================================
 * L5 Tests: Response Computation Accuracy
 * ========================================================================== */

static void test_first_order_step_response(void) {
    TEST("First-order step response accuracy");

    FirstOrderModel fo = {5.0, 2.0};

    /* At t=0 */
    assert(fabs(first_order_step(&fo, 0.0) - 0.0) < TOL);

    /* At t->inf */
    double y_large = first_order_step(&fo, 1e6);
    assert(fabs(y_large - 5.0) < 1e-6);

    /* Exact value at t=tau */
    double y_tau = first_order_step(&fo, 2.0);
    assert(fabs(y_tau - 5.0 * (1.0 - 1.0/M_E)) < TOL);

    /* Negative time */
    assert(fabs(first_order_step(&fo, -1.0)) < TOL);

    PASS();
}

static void test_second_order_impulse_energy(void) {
    TEST("Second-order impulse energy (Parseval)");

    /* For zeta=0.5, wn=2, K=1: E = K^2*wn/(4*zeta) = 1*2/(4*0.5) = 1 */
    SecondOrderModel so = {1.0, 0.5, 2.0};
    double E = second_order_impulse_energy(&so);
    assert(fabs(E - 1.0) < TOL);

    /* Critically damped: zeta=1, E = K^2*wn/4 */
    SecondOrderModel cd = {2.0, 1.0, 3.0};
    double E_cd = second_order_impulse_energy(&cd);
    assert(fabs(E_cd - 4.0*3.0/4.0) < TOL);  /* 4*3/4 = 3 */

    PASS();
}

static void test_step_response_overshoot(void) {
    TEST("Step response overshoot detection");

    SecondOrderModel ud = {1.0, 0.3, 1.0};  /* Underdamped */
    SecondOrderModel cd = {1.0, 1.0, 1.0};  /* Critically damped */

    assert(step_response_has_overshoot(&ud) == 1);
    assert(step_response_has_overshoot(&cd) == 0);

    PASS();
}

static void test_polynomial_operations(void) {
    TEST("Polynomial evaluation and division");

    /* P(s) = 1 + 2s + 3s^2; P(2) = 1 + 4 + 12 = 17 */
    double coeffs[3] = {1.0, 2.0, 3.0};
    assert(fabs(polynomial_eval_real(coeffs, 2, 2.0) - 17.0) < TOL);

    /* Derivative at zero: P'(0) = coeffs[1]*1! = 2 */
    assert(fabs(polynomial_derivative_at_zero(coeffs, 2, 1) - 2.0) < TOL);

    /* Polynomial division: s^2 / (s+1) = s-1 + 1/(s+1)
     * num=[0,0,1](deg2), den=[1,1](deg1) -> quotient=[-1,1], remainder=[1] */
    double num[3] = {0.0, 0.0, 1.0};
    double den[2] = {1.0, 1.0};
    double quot[2], rem[2];
    polynomial_division(num, 2, den, 1, quot, rem);
    assert(fabs(quot[0] - (-1.0)) < TOL);
    assert(fabs(quot[1] - 1.0) < TOL);

    PASS();
}

static void test_matrix_exponential(void) {
    TEST("Matrix exponential exp(A*t)");

    /* exp(0*t) = I */
    double A[4] = {0.0, 0.0, 0.0, 0.0};
    double expAt[4];
    matrix_exponential(A, 2, 1.0, expAt);
    assert(fabs(expAt[0] - 1.0) < TOL);
    assert(fabs(expAt[1] - 0.0) < TOL);
    assert(fabs(expAt[2] - 0.0) < TOL);
    assert(fabs(expAt[3] - 1.0) < TOL);

    /* exp(diag(1,2)*t) = diag(exp(t), exp(2t)) */
    double B[4] = {1.0, 0.0, 0.0, 2.0};
    matrix_exponential(B, 2, 1.0, expAt);
    assert(fabs(expAt[0] - exp(1.0)) < TOL);
    assert(fabs(expAt[3] - exp(2.0)) < TOL);

    PASS();
}

static void test_linear_solver(void) {
    TEST("Linear system solver");

    /* [[1,2],[3,4]] * x = [5,11] -> x = [1,2] */
    double A[4] = {1.0, 2.0, 3.0, 4.0};
    double b[2] = {5.0, 11.0};
    double x[2];
    int status = linear_solve(A, b, x, 2);
    assert(status == 0);
    assert(fabs(x[0] - 1.0) < TOL);
    assert(fabs(x[1] - 2.0) < TOL);

    PASS();
}

/* ==========================================================================
 * L6 Tests: Canonical Problem Verification
 * ========================================================================== */

static void test_second_order_underdamped_step(void) {
    TEST("Second-order underdamped step response");

    /* zeta=0.5, wn=1, K=2
     * Peak time: t_p = pi/(wn*sqrt(1-zeta^2)) = pi/(1*0.866) = 3.6276
     * Overshoot: M_p = exp(-pi*zeta/sqrt(1-zeta^2)) = exp(-pi*0.5/0.866) = 0.1630
     * So peak = K*(1+M_p) = 2*1.1630 = 2.3260
     */
    SecondOrderModel so = {2.0, 0.5, 1.0};

    double t_p_expected = M_PI / sqrt(0.75);
    double y_peak = second_order_step(&so, t_p_expected);

    double M_p = (y_peak - 2.0) / 2.0;
    double M_p_expected = exp(-M_PI * 0.5 / sqrt(0.75));

    assert(fabs(M_p - M_p_expected) < 1e-4);

    PASS();
}

static void test_second_order_critical_step(void) {
    TEST("Second-order critically damped step");

    /* zeta=1, wn=1, K=1: y(t) = 1 - (1+t)*exp(-t)
     * At t=1: y(1) = 1 - 2*exp(-1) = 0.26424
     * At t=5: y(5) = 1 - 6*exp(-5) = 0.95957
     */
    SecondOrderModel so = {1.0, 1.0, 1.0};

    double y1 = second_order_step(&so, 1.0);
    assert(fabs(y1 - (1.0 - 2.0 * exp(-1.0))) < TOL);

    double y5 = second_order_step(&so, 5.0);
    assert(fabs(y5 - (1.0 - 6.0 * exp(-5.0))) < TOL);

    PASS();
}

static void test_first_order_time_constant(void) {
    TEST("First-order time to fraction");

    FirstOrderModel fo = {1.0, 2.0};

    /* Time to reach 50% = -tau*ln(0.5) = -2*ln(0.5) = 1.38629 */
    double t50 = first_order_time_to_fraction(&fo, 0.5);
    assert(fabs(t50 - (-2.0 * log(0.5))) < TOL);

    /* Time to reach 90% = -tau*ln(0.1) = 4.60517 */
    double t90 = first_order_time_to_fraction(&fo, 0.9);
    assert(fabs(t90 - (-2.0 * log(0.1))) < TOL);

    PASS();
}

/* ==========================================================================
 * L7 Tests: Applications
 * ========================================================================== */

static void test_dc_motor_step(void) {
    TEST("DC motor step response");

    /* K=20, tau_e=0.001, tau_m=0.02: at t=0.2s (10*tau_m), near steady-state */
    double y = dc_motor_step_response(20.0, 0.001, 0.02, 0.2);
    assert(fabs(y - 20.0) < 0.01);  /* Should be very close after 10*tau_m */

    /* At t=0: should be 0 */
    assert(fabs(dc_motor_step_response(20.0, 0.001, 0.02, 0.0)) < TOL);

    PASS();
}

static void test_thermal_step(void) {
    TEST("Thermal step response");

    /* Q=10W, R_th=5K/W, C_th=100J/K -> tau=500s, K=50K steady-state rise
     * At t=tau: Delta_T = 10*5*(1-1/e) = 50*0.632 = 31.6 */
    double DT = thermal_step_response(10.0, 5.0, 100.0, 500.0);
    double expected = 10.0 * 5.0 * (1.0 - 1.0/M_E);
    assert(fabs(DT - expected) < 0.01);

    PASS();
}

static void test_vehicle_step(void) {
    TEST("Vehicle longitudinal step response");

    /* m=1300, b=50: K=0.02, tau=26
     * At t=26: v = 0.02*(1-1/e) ≈ 0.01265 m/s per Newton */
    double v = vehicle_longitudinal_step(1300.0, 50.0, 26.0);
    double expected = (1.0/50.0) * (1.0 - 1.0/M_E);
    assert(fabs(v - expected) < TOL);

    PASS();
}

/* ==========================================================================
 * L5 Tests: Convolution
 * ========================================================================== */

static void test_convolution_commutativity(void) {
    TEST("Convolution commutativity");

    int N = 100;
    double dt = 0.01;
    double *h = (double *)malloc((size_t)N * sizeof(double));
    double *u = (double *)malloc((size_t)N * sizeof(double));

    /* Simple exponential h and step u */
    for (int i = 0; i < N; i++) {
        double t = (double)i * dt;
        h[i] = exp(-t);
        u[i] = 1.0;  /* Step input */
    }

    int ok = convolution_commutative_check(h, u, N, dt, 1e-6);
    assert(ok == 1);

    free(h); free(u);
    PASS();
}

/* ==========================================================================
 * L5 Tests: Response Metrics
 * ========================================================================== */

static void test_response_metrics_second_order(void) {
    TEST("Response metrics for second-order system");

    SecondOrderModel so = {2.0, 0.4, 5.0};
    ResponseMetrics m;
    second_order_theoretical_metrics(&so, 0.02, &m);

    /* Peak time: t_p = pi/(wn*sqrt(1-zeta^2)) = pi/(5*0.9165) = 0.685 */
    double tp_expected = M_PI / (5.0 * sqrt(0.84));
    assert(fabs(m.peak_time - tp_expected) < 1e-3);

    /* Overshoot: exp(-pi*zeta/sqrt(1-zeta^2)) = exp(-pi*0.4/0.9165) = 0.254 */
    double ov_expected = exp(-M_PI * 0.4 / sqrt(0.84));
    assert(fabs(m.overshoot - ov_expected) < 1e-4);

    /* Steady-state = K = 2 */
    assert(fabs(m.steady_state - 2.0) < TOL);

    assert(validate_metrics(&m) == 1);

    PASS();
}

static void test_ziegler_nichols(void) {
    TEST("Ziegler-Nichols step response method");

    /* Generate a first-order step response and apply ZN method */
    FirstOrderModel fo = {3.0, 10.0};
    ResponseTrajectory *traj = first_order_step_trajectory(&fo, 50.0, 200);
    assert(traj != NULL);

    double K, L, T, Kp, Ti, Td;
    zieglier_nichols_step_method(traj, &K, &L, &T, &Kp, &Ti, &Td);

    /* For pure first-order (no dead time), L should be small */
    assert(K > 0.0);
    assert(T > 0.0);

    response_trajectory_free(traj);
    PASS();
}

/* ==========================================================================
 * L5 Tests: Discrete-Time Response
 * ========================================================================== */

static void test_discrete_step_response(void) {
    TEST("Discrete step response -- first-order");

    /* H(z) = 0.1*z^{-1} / (1 - 0.9*z^{-1})
     * This is the ZOH of G(s) = 1/(10s+1) with Ts = -(10)*ln(0.9) ≈ 1.0536 */
    DiscreteTransferFunction dtf;
    double num[2] = {0.0, 0.1};
    double den[2] = {1.0, -0.9};
    dtf.num = num; dtf.num_deg = 1;
    dtf.den = den; dtf.den_deg = 1;

    double s[50];
    discrete_step_response(&dtf, 50, s);

    /* Final value: H(1) = 0.1/(1-0.9) = 1.0 */
    double fv = discrete_final_value(&dtf);
    assert(fabs(fv - 1.0) < TOL);

    /* Check convergence */
    assert(fabs(s[49] - 1.0) < 0.01);

    PASS();
}

static void test_discrete_stability(void) {
    TEST("Discrete Jury stability test");

    /* H(z) = 1/(1 - 0.5*z^{-1}) -> pole at z=0.5, stable */
    DiscreteTransferFunction dtf_stable;
    double nums[1] = {1.0};
    double dens[2] = {1.0, -0.5};
    dtf_stable.num = nums; dtf_stable.num_deg = 0;
    dtf_stable.den = dens; dtf_stable.den_deg = 1;

    assert(discrete_is_stable(&dtf_stable) == 1);

    /* H(z) = 1/(1 - 1.5*z^{-1}) -> pole at z=1.5, unstable */
    DiscreteTransferFunction dtf_unstable;
    double denu[2] = {1.0, -1.5};
    dtf_unstable.num = nums; dtf_unstable.num_deg = 0;
    dtf_unstable.den = denu; dtf_unstable.den_deg = 1;

    assert(discrete_is_stable(&dtf_unstable) == 0);

    PASS();
}

/* ==========================================================================
 * L2 Tests: State-Space Conversion
 * ========================================================================== */

static void test_state_space_conversion(void) {
    TEST("Second-order to state-space conversion");

    SecondOrderModel so = {3.0, 0.7, 4.0};
    StateSpaceModel ss = {0};
    second_order_to_state_space(&so, &ss);

    /* Check dimensions */
    assert(ss.n_states == 2);
    assert(ss.n_inputs == 1);
    assert(ss.n_outputs == 1);

    /* Check A matrix: [[0, 1], [-wn^2, -2*zeta*wn]] = [[0,1], [-16, -5.6]] */
    assert(fabs(ss.A[0] - 0.0) < TOL);
    assert(fabs(ss.A[1] - 1.0) < TOL);
    assert(fabs(ss.A[2] - (-16.0)) < TOL);
    assert(fabs(ss.A[3] - (-5.6)) < TOL);

    /* Check B: [[0], [K*wn^2]] = [[0], [48]] */
    assert(fabs(ss.B[0]) < TOL);
    assert(fabs(ss.B[1] - 48.0) < TOL);

    /* Check C: [[1, 0]] */
    assert(fabs(ss.C[0] - 1.0) < TOL);
    assert(fabs(ss.C[1]) < TOL);

    state_space_free_matrices(&ss);
    PASS();
}

/* ==========================================================================
 * L4 Test: Impulse-Step Integral Relationship via Numerical
 * ========================================================================== */

static void test_impulse_integral_equals_step(void) {
    TEST("Integral of impulse response equals step response");

    FirstOrderModel fo = {2.0, 0.5};
    int N = 500;
    double t_final = 3.0;

    ResponseTrajectory *imp_traj = first_order_impulse_trajectory(&fo, t_final, N);
    assert(imp_traj != NULL);

    ResponseTrajectory *step_from_imp = step_from_impulse_integration(imp_traj);
    assert(step_from_imp != NULL);

    int mid_idx = N / 2;
    double t_mid = step_from_imp->data[mid_idx].t;
    double y_analytical = first_order_step(&fo, t_mid);
    double y_numerical = step_from_imp->data[mid_idx].y;

    assert(fabs(y_numerical - y_analytical) < 0.02);  /* 2% tolerance */

    response_trajectory_free(imp_traj);
    response_trajectory_free(step_from_imp);
    PASS();
}

/* ==========================================================================
 * Main test runner
 * ========================================================================== */

int main(void) {
    printf("\n=== mini-impulse-step-response: Test Suite ===\n\n");

    printf("--- L1: Core Definitions ---\n");
    test_damping_classification();
    test_second_order_poles();
    test_design_formulas();

    printf("\n--- L4: Mathematical Theorems ---\n");
    test_impulse_step_relationship();
    test_final_value_theorem();
    test_first_order_63_percent();
    test_stability_check();

    printf("\n--- L5: Response Computation ---\n");
    test_first_order_step_response();
    test_second_order_impulse_energy();
    test_step_response_overshoot();
    test_polynomial_operations();
    test_matrix_exponential();
    test_linear_solver();

    printf("\n--- L6: Canonical Problems ---\n");
    test_second_order_underdamped_step();
    test_second_order_critical_step();
    test_first_order_time_constant();

    printf("\n--- L7: Applications ---\n");
    test_dc_motor_step();
    test_thermal_step();
    test_vehicle_step();

    printf("\n--- L5: Convolution ---\n");
    test_convolution_commutativity();

    printf("\n--- L5: Response Metrics ---\n");
    test_response_metrics_second_order();
    test_ziegler_nichols();

    printf("\n--- L5: Discrete-Time Response ---\n");
    test_discrete_step_response();
    test_discrete_stability();

    printf("\n--- L2: State-Space Conversion ---\n");
    test_state_space_conversion();

    printf("\n--- L4: Impulse-Step Integral Relationship ---\n");
    test_impulse_integral_equals_step();

    printf("\n===========================================\n");
    printf("Tests run:    %d\n", tests_run);
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_run - tests_passed);
    printf("===========================================\n");

    return (tests_passed == tests_run) ? 0 : 1;
}
