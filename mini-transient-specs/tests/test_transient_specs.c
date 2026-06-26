/*============================================================================
 * test_transient_specs.c — Assert-based tests for transient specs library
 *
 * Tests core API functions: param extraction, spec computation,
 * stability (Routh-Hurwitz), IVT/FVT, design methods, performance
 * indices, sensitivity, and advanced topics.
 *============================================================================*/

#include "transient_specs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define TEST_EPS 1e-6
#define TEST_PI  3.14159265358979323846

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASSED\n"); } while(0)
#define FAIL(msg) do { printf("FAILED: %s\n", msg); } while(0)
#define CHECK(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)
#define CHECK_CLOSE(a, b, tol, msg) \
    do { if (fabs((a)-(b)) > (tol)) { \
        printf("FAILED: %s (%.6f vs %.6f)\n", msg, a, b); return; } \
    } while(0)

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * Test 1: Parameter Extraction
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

static void test_param_extraction(void)
{
    TEST("param_extraction_underdamped");
    second_order_params_t p;
    int ret = second_order_from_zeta_wn(0.5, 2.0, 1.0, &p);
    CHECK(ret == 0, "return code");
    CHECK(p.type == RESPONSE_UNDERDAMPED, "type");
    CHECK_CLOSE(p.damping_ratio, 0.5, TEST_EPS, "zeta");
    CHECK_CLOSE(p.natural_freq, 2.0, TEST_EPS, "wn");
    CHECK_CLOSE(p.damped_freq, 2.0*sqrt(0.75), TEST_EPS, "wd");
    CHECK_CLOSE(p.time_constant, 1.0, TEST_EPS, "tau");
    PASS();

    TEST("param_extraction_critically_damped");
    ret = second_order_from_zeta_wn(1.0, 3.0, 2.0, &p);
    CHECK(ret == 0, "return code");
    CHECK(p.type == RESPONSE_CRITICALLY_DAMPED, "type");
    CHECK_CLOSE(p.damped_freq, 0.0, TEST_EPS, "wd=0");
    PASS();

    TEST("param_extraction_overdamped");
    ret = second_order_from_zeta_wn(2.0, 2.0, 1.0, &p);
    CHECK(ret == 0, "return code");
    CHECK(p.type == RESPONSE_OVERDAMPED, "type");
    CHECK(p.damped_freq > 0.0, "wd>0 for overdamped");
    PASS();

    TEST("param_extraction_null");
    ret = second_order_from_zeta_wn(0.5, 2.0, 1.0, NULL);
    CHECK(ret == -1, "null params");
    PASS();
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * Test 2: Transient Spec Computation
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

static void test_spec_computation(void)
{
    TEST("specs_underdamped");
    second_order_params_t p;
    second_order_from_zeta_wn(0.5, 2.0, 1.0, &p);
    transient_specs_t s = compute_specs_second_order(&p);

    /* Expected: zeta=0.5, wn=2, wd=sqrt(3)=1.732
     * %OS = 100*exp(-pi*0.5/sqrt(0.75)) = 100*exp(-1.814) = 16.30%
     * tp = pi/1.732 = 1.814s
     * ts(2%) = 4/(0.5*2) = 4.0s
     * tr ~ 1.8/(0.5*2) = 1.8s */
    CHECK_CLOSE(s.percent_overshoot, 16.3034, 0.1, "%OS");
    CHECK_CLOSE(s.peak_time, TEST_PI/sqrt(3.0), 0.01, "tp");
    CHECK_CLOSE(s.settling_time_2pct, 4.0, 0.01, "ts2");
    CHECK_CLOSE(s.settling_time_5pct, 3.0, 0.01, "ts5");
    CHECK(s.steady_state_error == 0.0, "ess=0");
    PASS();

    TEST("specs_critically_damped");
    second_order_from_zeta_wn(1.0, 2.0, 1.0, &p);
    s = compute_specs_second_order(&p);
    CHECK(s.percent_overshoot == 0.0, "no overshoot");
    CHECK(!isfinite(s.peak_time), "no peak time");
    CHECK(s.settling_time_2pct > 0.0, "finite ts");
    PASS();

    TEST("specs_overdamped");
    second_order_from_zeta_wn(2.0, 2.0, 1.0, &p);
    s = compute_specs_second_order(&p);
    CHECK(s.percent_overshoot == 0.0, "no overshoot");
    CHECK(s.settling_time_2pct > 0.0, "finite ts");
    PASS();
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * Test 3: Routh-Hurwitz Stability
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

static void test_routh_hurwitz(void)
{
    TEST("routh_stable_3rd_order");
    double coeff_stable[] = {1.0, 6.0, 11.0, 6.0}; /* (s+1)(s+2)(s+3) */
    routh_table_t table;
    int ret = routh_hurwitz(coeff_stable, 3, &table);
    CHECK(ret == 0, "routh construction");
    CHECK(table.is_stable == 1, "stable");
    CHECK(table.sign_changes == 0, "no sign changes");
    PASS();

    TEST("routh_unstable");
    double coeff_unstable[] = {1.0, 1.0, -2.0, 1.0}; /* has RHP root */
    ret = routh_hurwitz(coeff_unstable, 3, &table);
    CHECK(ret == 0, "routh construction");
    CHECK(table.is_stable == 0, "unstable");
    CHECK(table.sign_changes > 0, "has sign changes");
    PASS();

    TEST("routh_null_input");
    ret = routh_hurwitz(NULL, 3, &table);
    CHECK(ret == -1, "null coeff");
    PASS();
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * Test 4: Final Value & Initial Value Theorems
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

static void test_value_theorems(void)
{
    TEST("final_value_theorem");
    /* G(s) = 2/(s+1), step => y(inf) = 2 */
    double num[] = {2.0};
    double den[] = {1.0, 1.0};
    double fv = final_value_theorem(num, 0, den, 1);
    CHECK_CLOSE(fv, 2.0, TEST_EPS, "FVT step");
    PASS();

    TEST("initial_value_theorem");
    /* G(s) = (s+2)/(s+3), proper => y(0+) = 1/1 = 1 */
    double num2[] = {2.0, 1.0};
    double den2[] = {3.0, 1.0};
    double iv = initial_value_theorem(num2, 1, den2, 1);
    CHECK_CLOSE(iv, 1.0, TEST_EPS, "IVT proper");
    PASS();

    TEST("initial_value_strictly_proper");
    /* G(s) = 2/(s+1), strictly proper => y(0+) = 0 */
    iv = initial_value_theorem(num, 0, den, 1);
    CHECK_CLOSE(iv, 0.0, TEST_EPS, "IVT strictly proper");
    PASS();
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * Test 5: Design Methods
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

static void test_design_methods(void)
{
    TEST("design_from_OS_ts");
    second_order_params_t p;
    int ret = design_from_OS_ts(16.3, 4.0, &p);
    CHECK(ret == 0, "design success");
    CHECK_CLOSE(p.damping_ratio, 0.5, 0.01, "zeta ~0.5");
    CHECK_CLOSE(p.natural_freq, 2.0, 0.01, "wn ~2");
    PASS();

    TEST("zeta_from_percent_overshoot");
    double z = zeta_from_percent_overshoot(16.3034);
    CHECK_CLOSE(z, 0.5, 0.01, "inverse OS");
    z = zeta_from_percent_overshoot(0.0);
    CHECK_CLOSE(z, 1.0, TEST_EPS, "0 OS -> zeta=1");
    z = zeta_from_percent_overshoot(100.0);
    CHECK_CLOSE(z, 0.0, TEST_EPS, "100 OS -> zeta=0");
    PASS();
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * Test 6: Dominant Pole Analysis
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

static void test_dominant_pole(void)
{
    TEST("dominant_pole_detection");
    pole_zero_descr_t pz;
    pole_zero_init(&pz, 4, 0);

    /* Poles: -1+j, -1-j (dominant pair), -5, -10 */
    pz.poles_real[0] = -1.0; pz.poles_imag[0] =  1.0;
    pz.poles_real[1] = -1.0; pz.poles_imag[1] = -1.0;
    pz.poles_real[2] = -5.0; pz.poles_imag[2] =  0.0;
    pz.poles_real[3] = -10.0; pz.poles_imag[3] = 0.0;
    pz.dc_gain = 1.0;

    dominant_pole_t dom;
    int ret = find_dominant_poles(&pz, &dom);
    CHECK(ret == 1, "found dominant");
    CHECK(dom.is_valid == 1, "dominance valid (ratio=5)");
    CHECK_CLOSE(dom.dominance_ratio, 5.0, 0.01, "dominance ratio");
    CHECK(dom.reduced_model.type == RESPONSE_UNDERDAMPED, "underdamped reduced");

    pole_zero_free(&pz);
    PASS();
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * Test 7: Performance Indices
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

static void test_performance_indices(void)
{
    TEST("performance_indices");
    second_order_params_t p;
    second_order_from_zeta_wn(0.5, 2.0, 1.0, &p);

    double iae = compute_IAE(&p);
    double ise = compute_ISE(&p);
    double itae = compute_ITAE(&p);
    double itse = compute_ITSE(&p);

    CHECK(iae > 0.0, "IAE > 0");
    CHECK(ise > 0.0, "ISE > 0");
    CHECK(itae > 0.0, "ITAE > 0");
    CHECK(itse > 0.0, "ITSE > 0");

    /* ISE should be smaller than IAE for same system (squared errors
     * are smaller when error < 1) */
    CHECK(ise < iae, "ISE < IAE for unit step");
    /* ITAE and IAE are both positive and finite */
    CHECK(itae > 0.0, "ITAE > 0");
    PASS();
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * Test 8: Sensitivity Analysis
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

static void test_sensitivity(void)
{
    TEST("sensitivity_analysis");
    second_order_params_t p;
    second_order_from_zeta_wn(0.5, 2.0, 1.0, &p);

    spec_sensitivity_t sens = compute_sensitivity(&p);
    CHECK(sens.d_ts_d_zeta < 0.0, "dts/dzeta negative");
    CHECK(sens.d_ts_d_wn < 0.0, "dts/dwn negative");
    CHECK(sens.d_tr_d_zeta < 0.0, "dtr/dzeta negative");

    /* Perturbation test: increase zeta by 0.1 => tr and ts decrease */
    transient_specs_t perturbed = specs_under_perturbation(&p, 0.1, 0.0);
    transient_specs_t nominal = compute_specs_second_order(&p);
    CHECK(perturbed.settling_time_2pct < nominal.settling_time_2pct,
          "ts decreases with higher zeta");
    PASS();
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * Test 9: Advanced Topics
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

static void test_advanced(void)
{
    TEST("delay_margin");
    second_order_params_t p;
    second_order_from_zeta_wn(0.707, 10.0, 1.0, &p);
    double dm = delay_margin_from_specs(&p, 60.0); /* 60 degree phase margin */
    CHECK(dm > 0.0, "delay margin positive");
    PASS();

    TEST("bandwidth");
    double bw = bandwidth_from_specs(&p);
    CHECK(bw > 0.0, "bandwidth positive");
    /* For zeta=0.707, bw ~ wn */
    CHECK_CLOSE(bw, 10.0, 2.0, "bw approx wn for zeta=0.707");
    PASS();

    TEST("time_varying_zeta");
    transient_specs_t s = specs_time_varying_zeta(0.5, 0.01, 2.0, 1.0);
    /* Should be similar to constant zeta specs */
    CHECK(s.settling_time_2pct > 0.0, "finite ts for time-varying");
    PASS();
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * Test 10: First-Order Systems
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

static void test_first_order(void)
{
    TEST("first_order_step");
    double y = first_order_step_response(0.0, 1.0, 1.0);
    CHECK_CLOSE(y, 0.0, TEST_EPS, "y(0)=0");
    y = first_order_step_response(1.0, 1.0, 1.0);
    CHECK_CLOSE(y, 1.0 - exp(-1.0), TEST_EPS, "y(tau)");
    y = first_order_step_response(10.0, 1.0, 1.0);
    CHECK_CLOSE(y, 1.0, 0.001, "y(inf) ~1");
    PASS();

    TEST("first_order_specs");
    transient_specs_t s = compute_specs_first_order(0.5, 1.0);
    CHECK_CLOSE(s.settling_time_2pct, 2.0, TEST_EPS, "ts=4tau");
    CHECK_CLOSE(s.rise_time, 0.5*log(9.0), TEST_EPS, "tr=tau*ln9");
    CHECK(s.percent_overshoot == 0.0, "no overshoot");
    CHECK(s.num_oscillations == 0, "no oscillations");
    PASS();
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * Test 11: Pole-Zero Utilities
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

static void test_pole_zero_utils(void)
{
    TEST("pole_zero_init_free");
    pole_zero_descr_t pz;
    int ret = pole_zero_init(&pz, 3, 2);
    CHECK(ret == 0, "init success");
    CHECK(pz.num_poles == 3, "num_poles");
    CHECK(pz.num_zeros == 2, "num_zeros");
    CHECK(pz.poles_real != NULL, "poles_real allocated");
    CHECK(pz.zeros_real != NULL, "zeros_real allocated");

    /* Write and read back */
    pz.poles_real[0] = -1.0; pz.poles_imag[0] = 2.0;
    pz.poles_real[1] = -1.0; pz.poles_imag[1] = -2.0;
    pz.poles_real[2] = -5.0; pz.poles_imag[2] = 0.0;
    pz.zeros_real[0] = -2.0; pz.zeros_imag[0] = 0.0;
    pz.zeros_real[1] = -3.0; pz.zeros_imag[1] = 1.0;

    CHECK_CLOSE(pz.poles_real[0], -1.0, TEST_EPS, "pole 0 real");
    CHECK_CLOSE(pz.poles_imag[0], 2.0, TEST_EPS, "pole 0 imag");

    pole_zero_free(&pz);
    CHECK(pz.poles_real == NULL, "freed");
    CHECK(pz.num_poles == 0, "zeroed");
    PASS();
}

/*━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 * Main test runner
 *━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━*/

int main(void)
{
    printf("=== transient_specs Test Suite ===\n\n");

    test_param_extraction();
    test_spec_computation();
    test_routh_hurwitz();
    test_value_theorems();
    test_design_methods();
    test_dominant_pole();
    test_performance_indices();
    test_sensitivity();
    test_advanced();
    test_first_order();
    test_pole_zero_utils();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
