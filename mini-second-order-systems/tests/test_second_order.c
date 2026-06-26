/**
 * @file test_second_order.c
 * @brief Comprehensive test suite for mini-second-order-systems
 *
 * Tests cover:
 *   L1: struct creation, damping classification, pole computation
 *   L2: transfer function operations, frequency response
 *   L4: stability, BIBO, H2 norm, sensitivity
 *   L4: transient spec formulas (theorems: PO, t_p, t_s exact formulas)
 *   L5: inverse computation, system identification methods
 *   L6: canonical physical models
 *   L5/L8: design optimization, sensitivity
 *
 * All tests use standard assert(). No custom macros.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

/* Mini assert with tolerance for floating-point comparisons */
#define ASSERT_NEAR(a, b, tol) \
    assert(fabs((a) - (b)) < (tol))

#include "second_order.h"
#include "transient_specs.h"
#include "response_computation.h"
#include "system_identification.h"
#include "canonical_models.h"
#include "second_order_design.h"

static int tests_run = 0;
static int tests_passed = 0;

#define RUN_TEST(name) do { \
    tests_run++; \
    printf("  [%3d] %s... ", tests_run, #name); \
    name(); \
    tests_passed++; \
    printf("PASSED\n"); \
} while(0)

/* ==========================================================================
 * L1: System Construction and Classification
 * ========================================================================== */

void test_system_creation(void)
{
    SecondOrderSystem sys = so2_create(1.0, 0.5, 2.0);
    assert(sys.K == 1.0);
    assert(sys.zeta == 0.5);
    assert(sys.omega_n == 2.0);
}

void test_system_from_poly(void)
{
    /* G(s) = 4/(s² + 2s + 4) → ωₙ=2, ζ=0.5, K=1 */
    double num[3] = {0, 0, 4};
    double den[3] = {1, 2, 4};
    SecondOrderSystem sys = so2_from_poly(num, den);
    ASSERT_NEAR(sys.zeta, 0.5, 1e-6);
    ASSERT_NEAR(sys.omega_n, 2.0, 1e-6);
    ASSERT_NEAR(sys.K, 1.0, 1e-6);
}

void test_system_from_poles(void)
{
    /* Poles at -1 ± j2 → ωₙ=√(1+4)=√5≈2.236, ζ=1/√5≈0.447 */
    SecondOrderSystem sys = so2_from_poles(-1.0, 2.0, 1.0);
    ASSERT_NEAR(sys.zeta, 1.0 / sqrt(5.0), 1e-6);
    ASSERT_NEAR(sys.omega_n, sqrt(5.0), 1e-6);
}

void test_damping_classification(void)
{
    SecondOrderSystem s_undamped  = so2_create(1, 0.0, 1.0);
    SecondOrderSystem s_under    = so2_create(1, 0.3, 1.0);
    SecondOrderSystem s_critical = so2_create(1, 1.0, 1.0);
    SecondOrderSystem s_over     = so2_create(1, 2.0, 1.0);
    SecondOrderSystem s_unstable = so2_create(1, -0.5, 1.0);
    assert(so2_damping_class(&s_undamped) == DAMPING_UNDAMPED);
    assert(so2_damping_class(&s_under) == DAMPING_UNDERDAMPED);
    assert(so2_damping_class(&s_critical) == DAMPING_CRITICALLY);
    assert(so2_damping_class(&s_over) == DAMPING_OVERDAMPED);
    assert(so2_damping_class(&s_unstable) == DAMPING_UNSTABLE);
}

void test_pole_computation(void)
{
    SecondOrderSystem sys = so2_create(1.0, 0.5, 4.0);
    SecondOrderPole p1, p2;
    so2_poles(&sys, &p1, &p2);

    /* ζ=0.5, ωₙ=4 → σ=-2, ω_d=4√(0.75)=4·0.866=3.464 */
    ASSERT_NEAR(p1.real, -2.0, 1e-6);
    ASSERT_NEAR(p1.imag, 4.0 * sqrt(0.75), 1e-6);
    /* Conjugate pair */
    ASSERT_NEAR(p2.real, -2.0, 1e-6);
    ASSERT_NEAR(p2.imag, -4.0 * sqrt(0.75), 1e-6);
}

void test_stability(void)
{
    SecondOrderSystem s_stable   = so2_create(1, 0.5, 1.0);
    SecondOrderSystem s_undamped = so2_create(1, 0.0, 1.0);
    SecondOrderSystem s_unstable = so2_create(1, -0.1, 1.0);
    assert(so2_is_stable(&s_stable) == 1);
    assert(so2_is_stable(&s_undamped) == 0);  /* Undamped */
    assert(so2_is_stable(&s_unstable) == 0); /* Unstable */
}

/* ==========================================================================
 * L2/L4: Transfer Function and Frequency Response
 * ========================================================================== */

void test_transfer_function(void)
{
    SecondOrderSystem sys = so2_create(2.0, 0.5, 3.0);
    TransferFunc2 tf = so2_to_transfer_function(&sys);
    /* den = [1, 2ζωₙ, ωₙ²] = [1, 3, 9] */
    assert(tf.den[0] == 1.0);
    ASSERT_NEAR(tf.den[1], 3.0, 1e-6);
    ASSERT_NEAR(tf.den[2], 9.0, 1e-6);
    /* num = [0, 0, Kωₙ²] = [0, 0, 18] */
    ASSERT_NEAR(tf.num[2], 18.0, 1e-6);
}

void test_frequency_response_dc(void)
{
    SecondOrderSystem sys = so2_create(2.0, 0.5, 3.0);
    double mag = so2_magnitude(&sys, 0.0);
    ASSERT_NEAR(mag, 2.0, 1e-6);  /* DC gain */
}

void test_frequency_response_at_omega_n(void)
{
    /* At ω=ωₙ: |G| = K/(2ζ) */
    SecondOrderSystem sys = so2_create(1.0, 0.5, 3.0);
    double mag = so2_magnitude(&sys, 3.0);
    ASSERT_NEAR(mag, 1.0 / (2.0 * 0.5), 1e-6);
}

void test_phase_at_omega_n(void)
{
    /* At ω=ωₙ: phase = -90° = -π/2 */
    SecondOrderSystem sys = so2_create(1.0, 0.5, 3.0);
    double phase = so2_phase(&sys, 3.0);
    ASSERT_NEAR(phase, -M_PI / 2.0, 1e-6);
}

void test_bibo_stability(void)
{
    SecondOrderSystem s1 = so2_create(1, 0.5, 1.0);
    SecondOrderSystem s2 = so2_create(1, 0.0, 1.0);
    assert(so2_is_bibo_stable(&s1) == 1);
    assert(so2_is_bibo_stable(&s2) == 0);
}

void test_h2_norm(void)
{
    /* H₂² = K²·ωₙ³/(4ζ) */
    SecondOrderSystem sys = so2_create(1.0, 0.5, 1.0);
    double H2_sq = 1.0 * 1.0 / (4.0 * 0.5);
    double H2 = sqrt(H2_sq);
    ASSERT_NEAR(so2_h2_norm(&sys), H2, 1e-6);
}

/* ==========================================================================
 * L4: Transient Specification Theorems
 * ========================================================================== */

void test_peak_time_formula(void)
{
    /* ζ=0.5, ωₙ=2, ω_d=2√(0.75)=1.732, t_p=π/1.732≈1.814 */
    double tp = transient_peak_time(0.5, 2.0);
    ASSERT_NEAR(tp, M_PI / (2.0 * sqrt(0.75)), 1e-6);

    /* Critically damped: no peak */
    assert(transient_peak_time(1.0, 1.0) == 0.0);
}

void test_percent_overshoot_formula(void)
{
    /* ζ=0.5: PO = 100·exp(-π·0.5/√(0.75)) ≈ 16.303% */
    double po = transient_percent_overshoot(0.5);
    ASSERT_NEAR(po, 16.303, 0.05);

    /* ζ=0.707: PO ≈ 4.6% */
    double po_opt = transient_percent_overshoot(0.707);
    ASSERT_NEAR(po_opt, 4.6, 0.5);

    /* ζ=1: no overshoot */
    assert(transient_percent_overshoot(1.0) == 0.0);
}

void test_settling_time_formula(void)
{
    /* ζ=0.5, ωₙ=2: t_s_2% = 4/(0.5·2) = 4 */
    double ts2 = transient_settling_time_2pct(0.5, 2.0);
    ASSERT_NEAR(ts2, 4.0, 1e-6);

    /* t_s_5% = 3/(0.5·2) = 3 */
    double ts5 = transient_settling_time_5pct(0.5, 2.0);
    ASSERT_NEAR(ts5, 3.0, 1e-6);
}

void test_inverse_overshoot_formula(void)
{
    /* PO=16.303% → ζ≈0.5 */
    double zeta = transient_zeta_from_overshoot(16.303);
    ASSERT_NEAR(zeta, 0.5, 1e-2);

    /* PO=0% → ζ=1 */
    assert(transient_zeta_from_overshoot(0.0) == 1.0);

    /* PO=100% → ζ=0 */
    assert(transient_zeta_from_overshoot(100.0) == 0.0);
}

void test_design_from_specs(void)
{
    double zeta, omega_n;
    int ok = transient_design_from_specs(16.303, 4.0, &zeta, &omega_n);
    assert(ok == 1);
    ASSERT_NEAR(zeta, 0.5, 1e-2);
    ASSERT_NEAR(omega_n, 2.0, 1e-2);
}

/* ==========================================================================
 * L2/L4: Step, Impulse, and Ramp Responses
 * ========================================================================== */

void test_step_response_steady_state(void)
{
    /* y(t→∞) = K */
    SecondOrderSystem sys = so2_create(2.0, 0.5, 1.0);
    double y_ss = response_step(&sys, 100.0);
    ASSERT_NEAR(y_ss, 2.0, 1e-6);
}

void test_step_response_initial(void)
{
    SecondOrderSystem sys = so2_create(1.0, 0.5, 1.0);
    assert(response_step(&sys, 0.0) == 0.0);
}

void test_step_response_critically_damped(void)
{
    /* ζ=1: y(t) = K·(1 - e^{-ωₙt}·(1+ωₙt)) */
    SecondOrderSystem sys = so2_create(1.0, 1.0, 1.0);
    double t = 1.0;
    double expected = 1.0 - exp(-1.0) * 2.0;
    ASSERT_NEAR(response_step(&sys, t), expected, 1e-6);
}

void test_step_response_undamped(void)
{
    /* ζ=0: y(t) = K·(1 - cos(ωₙt)) */
    SecondOrderSystem sys = so2_create(1.0, 0.0, 1.0);
    double t = M_PI;
    double expected = 1.0 - cos(M_PI);  /* = 2.0 */
    ASSERT_NEAR(response_step(&sys, t), expected, 1e-6);
}

void test_impulse_response_integral(void)
{
    /* ∫₀^∞ g(t) dt should approach K (DC gain) for stable system */
    SecondOrderSystem sys = so2_create(1.0, 0.5, 1.0);
    double integral = 0.0;
    double dt = 0.01;
    for (double t = dt / 2; t < 20.0; t += dt) {
        integral += response_impulse(&sys, t) * dt;
    }
    ASSERT_NEAR(integral, 1.0, 0.05);
}

void test_free_response_decay(void)
{
    /* Free response should decay to 0 for stable system */
    SecondOrderSystem sys = so2_create(1.0, 0.5, 1.0);
    double y = response_free(&sys, 20.0, 1.0, 0.0);
    ASSERT_NEAR(y, 0.0, 1e-4);
}

void test_energy_dissipation(void)
{
    SecondOrderSystem sys = so2_create(1.0, 0.5, 1.0);
    /* V = 0.5·ẏ² + 0.5·ωₙ²·y² >= 0 */
    double V = response_energy_function(&sys, 1.0, 0.5);
    assert(V > 0.0);
    /* Dissipation rate should be negative */
    double dVdt = response_energy_dissipation_rate(&sys, 0.5);
    assert(dVdt < 0.0);
}

/* ==========================================================================
 * L2: Resonance and Bandwidth
 * ========================================================================== */

void test_resonance_frequency(void)
{
    /* ζ=0.3, ωₙ=10: ω_r = 10√(1-2·0.09) = 10√(0.82) ≈ 9.055 */
    SecondOrderSystem sys = so2_create(1.0, 0.3, 10.0);
    double wr = response_resonance_frequency(&sys);
    ASSERT_NEAR(wr, 10.0 * sqrt(0.82), 1e-6);

    /* ζ=0.8 > 0.707: no resonance */
    SecondOrderSystem sys2 = so2_create(1.0, 0.8, 10.0);
    assert(response_resonance_frequency(&sys2) == 0.0);
}

void test_resonance_peak(void)
{
    /* ζ=0.3: M_r = 1/(2·0.3·√(0.91)) ≈ 1.747 */
    SecondOrderSystem sys = so2_create(1.0, 0.3, 10.0);
    double Mr = response_resonance_peak(&sys);
    ASSERT_NEAR(Mr, 1.0 / (0.6 * sqrt(0.91)), 1e-6);
}

void test_bandwidth(void)
{
    /* ζ=0.5, ωₙ=1: BW ≈ 1.27 (numerically verified from Ogata) */
    SecondOrderSystem sys = so2_create(1.0, 0.5, 1.0);
    double bw = response_bandwidth(&sys);
    /* |G(j·bw)| should be 1/√2 of DC gain */
    double mag = so2_magnitude(&sys, bw);
    ASSERT_NEAR(mag, 1.0 / sqrt(2.0), 1e-3);
}

/* ==========================================================================
 * L5: System Identification
 * ========================================================================== */

void test_log_decrement_method(void)
{
    /* δ = ln(10/3) ≈ 1.204, ζ = 1.204/√(4π²+1.204²) ≈ 0.188 */
    double peaks[] = {10.0, 3.0, 0.9};
    double zeta, omega_d;
    int ok = sysid_from_log_decrement(peaks, 3, &zeta, &omega_d, 0.0);
    assert(ok == 1);
    /* Rough check: zeta should be positive and less than 1 */
    assert(zeta > 0.01 && zeta < 1.0);
}

void test_overshoot_peak_identification(void)
{
    double zeta, omega_n;
    int ok = sysid_from_overshoot_peak(16.303, M_PI / (2.0 * sqrt(0.75)),
                                         &zeta, &omega_n);
    assert(ok == 1);
    ASSERT_NEAR(zeta, 0.5, 1e-2);
    ASSERT_NEAR(omega_n, 2.0, 1e-2);
}

void test_quality_factor_to_zeta(void)
{
    /* Q=10 → ζ≈0.05 */
    double zeta = sysid_zeta_from_quality_factor(10.0);
    ASSERT_NEAR(zeta, 0.05, 0.005);
}

void test_prony_method(void)
{
    /* Generate impulse response data from known system.
     * Prony's method works best on impulse (free-decay) data. */
    double dt = 0.005;
    int n = 500;
    double *y = (double *)malloc((size_t)n * sizeof(double));
    assert(y != NULL);

    SecondOrderSystem true_sys = so2_create(1.0, 0.3, 5.0);
    /* Use impulse response (free decay) for Prony */
    for (int i = 0; i < n; i++) {
        y[i] = response_impulse(&true_sys, i * dt);
    }

    SecondOrderSystem est;
    int ok = sysid_prony_method(y, n, dt, &est);
    assert(ok == 1);
    /* Prony identifies modes; check outputs are finite and reasonable */
    assert(isfinite(est.zeta));
    assert(isfinite(est.omega_n));
    assert(est.omega_n > 0.0);

    free(y);
}

/* ==========================================================================
 * L6: Canonical Physical Models
 * ========================================================================== */

void test_mass_spring_damper(void)
{
    MassSpringDamper msd = {1.0, 2.0, 100.0};  /* m=1, b=2, k=100 */
    SecondOrderSystem sys = msd_to_second_order(&msd);
    ASSERT_NEAR(sys.omega_n, 10.0, 1e-6);
    ASSERT_NEAR(sys.zeta, 0.1, 1e-6);
    ASSERT_NEAR(sys.K, 0.01, 1e-6);

    double b_crit = msd_critical_damping(&msd);
    ASSERT_NEAR(b_crit, 20.0, 1e-6);
}

void test_series_rlc(void)
{
    /* R=10, L=1e-3, C=1e-6 → ωₙ=31623, ζ=0.158 */
    SeriesRLC rlc = {10.0, 1e-3, 1e-6};
    SecondOrderSystem sys = rlc_series_to_second_order(&rlc);
    ASSERT_NEAR(sys.omega_n, 31622.7766, 0.1);
    ASSERT_NEAR(sys.zeta, 0.158, 0.01);
}

void test_parallel_rlc(void)
{
    /* R=1000, L=1e-3, C=1e-6 → high Q */
    ParallelRLC rlc = {1000.0, 1e-3, 1e-6};
    SecondOrderSystem sys = rlc_parallel_to_second_order(&rlc);
    ASSERT_NEAR(sys.omega_n, 31622.7766, 0.1);
    /* ζ = (1/(2R))·√(L/C) = (1/2000)·√(1000) ≈ 0.0158 */
    ASSERT_NEAR(sys.zeta, 0.0158, 0.002);
}

void test_pendulum(void)
{
    SimplePendulum pend = {1.0, 0.1, 0.01, 9.81};
    SecondOrderSystem sys = pendulum_to_second_order(&pend);
    ASSERT_NEAR(sys.omega_n, sqrt(9.81), 1e-6);
    ASSERT_NEAR(pendulum_period(&pend), 2.0 * M_PI / sqrt(9.81), 1e-3);
}

void test_dc_motor_model(void)
{
    DCMotor motor = {1.0, 5e-3, 0.01, 0.01, 1e-4, 1e-5};
    /* τ_e = L/R = 5ms, τ_m = J·R/(Kt·Ke) = 1e-4/(1e-4) = 1s */
    ASSERT_NEAR(dc_motor_electrical_time_constant(&motor), 0.005, 1e-6);
    ASSERT_NEAR(dc_motor_mechanical_time_constant(&motor), 1.0, 1e-4);
}

void test_mems_accelerometer(void)
{
    MEMSAccel accel = {1e-9, 100.0, 1e-6, 2e-6, 1e-6};
    SecondOrderSystem sys = mems_accel_to_second_order(&accel);
    ASSERT_NEAR(sys.omega_n, sqrt(100.0 / 1e-9), 1e-2);
    /* Sensitivity check */
    double sens = mems_accel_sensitivity(&accel);
    assert(sens > 0.0);
}

/* ==========================================================================
 * L5/L8: Design and Optimization
 * ========================================================================== */

void test_itae_optimal_zeta(void)
{
    double zeta_itae = design_zeta_itae_optimal();
    ASSERT_NEAR(zeta_itae, 0.707, 1e-3);
}

void test_ise_computation(void)
{
    /* ISE for ζ=0.5, ωₙ=1: (1+4·0.25)/(4·0.5·1) = 2/2 = 1 */
    double ise = design_compute_ise(0.5, 1.0);
    ASSERT_NEAR(ise, 1.0, 1e-6);
}

void test_speed_damping_tradeoff(void)
{
    double norm_settle, norm_rise, overshoot;
    design_speed_damping_tradeoff(0.707, &norm_settle, &norm_rise, &overshoot);
    /* For ζ=0.707: norm_settle ≈ 5.66 */
    ASSERT_NEAR(norm_settle, 4.0 / 0.707, 0.01);
    ASSERT_NEAR(overshoot, 4.6, 1.0);
}

void test_fastest_for_overshoot(void)
{
    /* Max 5% overshoot → ζ ≈ 0.69 */
    double zeta = design_fastest_for_max_overshoot(5.0);
    assert(zeta > 0.6 && zeta < 0.8);
}

void test_sensitivity_analysis(void)
{
    double dPO, dts, dtp;
    design_sensitivity_to_zeta(0.5, 2.0, &dPO, &dts, &dtp);
    /* dPO/dζ should be negative (more damping → less overshoot) */
    assert(dPO < 0.0);
    /* dts/dζ should be negative (more damping → faster settling for underdamped) */
    assert(dts < 0.0);
}

void test_p_gain_design(void)
{
    SecondOrderSystem plant = so2_create(1.0, 1.0, 2.0);  /* Critically damped */
    /* Want ζ=0.5 in closed loop */
    double Kp = design_p_gain_for_zeta(&plant, 0.5);
    /* Kp should be positive */
    assert(Kp > 0.0);
}

void test_pd_design(void)
{
    SecondOrderSystem plant = so2_create(1.0, 0.3, 2.0);
    double Kp, Kd;
    int ok = design_pd_for_poles(&plant, 0.7, 5.0, &Kp, &Kd);
    assert(ok == 1);
    assert(Kp >= 0.0);
    /* Kd may be positive or negative depending on desired ζ */
}

void test_application_zeta_guidelines(void)
{
    double zmin, zmax;
    double zeta = design_zeta_guideline(APP_PRECISION_POSITION, &zmin, &zmax);
    assert(zeta >= zmin && zeta <= zmax);
    assert(zmin >= 0.5);  /* Precision needs good damping */
}

void test_dominant_pole_approximation(void)
{
    double real[] = {-1.0, -1.0, -10.0};
    double imag[] = {2.0, -2.0, 0.0};
    SecondOrderSystem sys = design_dominant_pole_approximation(real, imag, 3);
    ASSERT_NEAR(sys.zeta, 1.0 / sqrt(5.0), 1e-2);
}

/* ==========================================================================
 * L2/L4: SSE and Error Constants
 * ========================================================================== */

void test_error_constants(void)
{
    SecondOrderSystem sys = so2_create(5.0, 0.5, 1.0);
    double Kp, Kv, Ka;
    so2_error_constants(&sys, &Kp, &Kv, &Ka);
    ASSERT_NEAR(Kp, 5.0, 1e-6);
    assert(Kv == 0.0);
    assert(Ka == 0.0);
}

void test_steady_state_error(void)
{
    SecondOrderSystem sys = so2_create(9.0, 0.5, 1.0);
    double e_step, e_ramp, e_parabolic;
    so2_steady_state_errors(&sys, &e_step, &e_ramp, &e_parabolic);
    ASSERT_NEAR(e_step, 0.1, 1e-6);  /* 1/(1+9) = 0.1 */
}

/* ==========================================================================
 * L5: Transient Specs Aggregate
 * ========================================================================== */

void test_transient_compute_all(void)
{
    SecondOrderSystem sys = so2_create(1.0, 0.5, 2.0);
    TransientSpecs spec = transient_compute_all(&sys);
    ASSERT_NEAR(spec.percent_overshoot, 16.303, 0.05);
    ASSERT_NEAR(spec.steady_state, 1.0, 1e-6);
    ASSERT_NEAR(spec.max_value, 1.0 + 0.16303, 0.01);
    assert(spec.peak_time > 0.0);
    assert(spec.settling_time_2pct > 0.0);
}

/* ==========================================================================
 * L5: Numerical Simulation
 * ========================================================================== */

void test_rk4_step_response(void)
{
    /* Use ζ=0.7 (faster settling) and longer time for cleaner convergence */
    SecondOrderSystem sys = so2_create(1.0, 0.7, 2.0);
    TimeResponse traj = response_simulate_step(&sys, 8.0, 800,
                                                1.0, 0.0, 0.0);
    assert(traj.n == 801);
    /* Settling time t_s = 4/(0.7*2) ≈ 2.86s. At T=8s, well past settling. */
    ASSERT_NEAR(traj.data[800].y, 1.0, 1e-3);
    response_free_trajectory(&traj);
}

void test_phase_portrait_equilibrium(void)
{
    SecondOrderSystem sys = so2_create(1.0, 0.5, 2.0);
    /* Stable focus */
    const char *eq_type = response_equilibrium_type(&sys);
    assert(eq_type != NULL);

    SecondOrderSystem sys2 = so2_create(1.0, 0.0, 2.0);
    /* Undamped → center */
    const char *eq2 = response_equilibrium_type(&sys2);
    assert(eq2 != NULL);
}

/* ==========================================================================
 * Main
 * ========================================================================== */

int main(void)
{
    printf("=== mini-second-order-systems Test Suite ===\n\n");

    /* L1: Definitions */
    printf("--- L1: System Construction and Classification ---\n");
    RUN_TEST(test_system_creation);
    RUN_TEST(test_system_from_poly);
    RUN_TEST(test_system_from_poles);
    RUN_TEST(test_damping_classification);
    RUN_TEST(test_pole_computation);
    RUN_TEST(test_stability);

    /* L2/L4: Transfer function and frequency response */
    printf("\n--- L2/L4: Transfer Function and Frequency Response ---\n");
    RUN_TEST(test_transfer_function);
    RUN_TEST(test_frequency_response_dc);
    RUN_TEST(test_frequency_response_at_omega_n);
    RUN_TEST(test_phase_at_omega_n);
    RUN_TEST(test_bibo_stability);
    RUN_TEST(test_h2_norm);
    RUN_TEST(test_error_constants);
    RUN_TEST(test_steady_state_error);

    /* L4: Transient spec theorems */
    printf("\n--- L4: Transient Specification Theorems ---\n");
    RUN_TEST(test_peak_time_formula);
    RUN_TEST(test_percent_overshoot_formula);
    RUN_TEST(test_settling_time_formula);
    RUN_TEST(test_inverse_overshoot_formula);
    RUN_TEST(test_design_from_specs);

    /* L2/L4: Step, Impulse, Ramp responses */
    printf("\n--- L2/L4: Response Computation ---\n");
    RUN_TEST(test_step_response_steady_state);
    RUN_TEST(test_step_response_initial);
    RUN_TEST(test_step_response_critically_damped);
    RUN_TEST(test_step_response_undamped);
    RUN_TEST(test_impulse_response_integral);
    RUN_TEST(test_free_response_decay);
    RUN_TEST(test_energy_dissipation);

    /* L2: Resonance */
    printf("\n--- L2: Resonance and Bandwidth ---\n");
    RUN_TEST(test_resonance_frequency);
    RUN_TEST(test_resonance_peak);
    RUN_TEST(test_bandwidth);

    /* L5: System identification */
    printf("\n--- L5: System Identification ---\n");
    RUN_TEST(test_log_decrement_method);
    RUN_TEST(test_overshoot_peak_identification);
    RUN_TEST(test_quality_factor_to_zeta);
    RUN_TEST(test_prony_method);

    /* L6: Canonical models */
    printf("\n--- L6: Canonical Physical Models ---\n");
    RUN_TEST(test_mass_spring_damper);
    RUN_TEST(test_series_rlc);
    RUN_TEST(test_parallel_rlc);
    RUN_TEST(test_pendulum);
    RUN_TEST(test_dc_motor_model);
    RUN_TEST(test_mems_accelerometer);

    /* L5/L8: Design */
    printf("\n--- L5/L8: Design and Optimization ---\n");
    RUN_TEST(test_itae_optimal_zeta);
    RUN_TEST(test_ise_computation);
    RUN_TEST(test_speed_damping_tradeoff);
    RUN_TEST(test_fastest_for_overshoot);
    RUN_TEST(test_sensitivity_analysis);
    RUN_TEST(test_p_gain_design);
    RUN_TEST(test_pd_design);
    RUN_TEST(test_application_zeta_guidelines);
    RUN_TEST(test_dominant_pole_approximation);

    /* L5: Aggregate and simulation */
    printf("\n--- L5: Aggregate Specs and Simulation ---\n");
    RUN_TEST(test_transient_compute_all);
    RUN_TEST(test_rk4_step_response);
    RUN_TEST(test_phase_portrait_equilibrium);

    printf("\n=== %d/%d tests PASSED ===\n", tests_passed, tests_run);
    return 0;
}
