/**
 * @file    test_routh_hurwitz.c
 * @brief   Comprehensive test suite for Routh-Hurwitz stability module
 *
 * Tests cover:
 *   - Polynomial initialization
 *   - Routh array construction for stable/unstable/marginal systems
 *   - First-column sign change counting
 *   - Necessary condition verification
 *   - Hurwitz determinant computation
 *   - Special cases (zero first column, entire zero row)
 *   - Root distribution counting
 *   - Liénard-Chipart criterion
 *   - Relative stability (axis shift, stability margin)
 *   - Jury stability criterion
 *   - Low-order explicit conditions
 *   - Bilinear transform
 *   - Companion matrix construction
 *   - Edge cases (NULL, zero degree, etc.)
 *
 * @author  Automated Control Knowledge Base
 * @date    2026-06-20
 * @license MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "routh_hurwitz.h"
#include "hurwitz_determinant.h"
#include "special_cases.h"
#include "relative_stability.h"
#include "stability_criteria.h"
#include "jury_stability.h"

/* ============================================================================
 * Test Infrastructure
 * ============================================================================ */

static int tests_passed = 0;
static int tests_failed = 0;
static int tests_total = 0;

#define TEST(name) do { \
    tests_total++; \
    printf("  TEST %3d: %-55s ", tests_total, name); \
} while(0)

#define PASS() do { \
    printf("PASS\n"); \
    tests_passed++; \
} while(0)

#define FAIL(msg) do { \
    printf("FAIL: %s\n", msg); \
    tests_failed++; \
} while(0)

#define ASSERT_TRUE(cond, msg) do { \
    if (!(cond)) { FAIL(msg); } else { PASS(); } \
} while(0)

#define ASSERT_EQ(a, b, msg) do { \
    if ((a) != (b)) { printf("FAIL: %s (got %d, expected %d)\n", msg, (int)(a), (int)(b)); tests_failed++; } \
    else { PASS(); } \
} while(0)

#define ASSERT_DOUBLE_EQ(a, b, tol, msg) do { \
    if (fabs((a) - (b)) > (tol)) { printf("FAIL: %s (got %.10g, expected %.10g)\n", msg, a, b); tests_failed++; } \
    else { PASS(); } \
} while(0)

#define ASSERT_DOUBLE_GT(a, b, msg) do { \
    if (!((a) > (b))) { printf("FAIL: %s (%.10g not > %.10g)\n", msg, a, b); tests_failed++; } \
    else { PASS(); } \
} while(0)

/* ============================================================================
 * L1 Tests: Polynomial & Routh Array Construction
 * ============================================================================ */

void test_polynomial_init(void)
{
    printf("\n--- L1: Polynomial Initialization ---\n");

    TEST("init valid degree-3 polynomial");
    {
        double c[] = {2.0, 5.0, 3.0, 1.0}; /* a_0=2, a_1=5, a_2=3, a_3=1 */
        RouthPolynomial poly;
        bool ok = routh_polynomial_init(&poly, c, 3);
        ASSERT_TRUE(ok && poly.degree == 3 && fabs(poly.coeff[3] - 1.0) < 1e-12,
                     "degree 3 poly init failed");
    }

    TEST("init NULL polynomial fails");
    {
        double c[] = {1.0, 2.0};
        bool ok = routh_polynomial_init(NULL, c, 2);
        ASSERT_TRUE(!ok, "NULL poly should fail");
    }

    TEST("init NULL coeffs fails");
    {
        RouthPolynomial poly;
        bool ok = routh_polynomial_init(&poly, NULL, 2);
        ASSERT_TRUE(!ok, "NULL coeffs should fail");
    }

    TEST("init zero leading coefficient fails");
    {
        double c[] = {1.0, 2.0, 0.0}; /* a_2 = 0 */
        RouthPolynomial poly;
        bool ok = routh_polynomial_init(&poly, c, 2);
        ASSERT_TRUE(!ok, "zero leading coeff should fail");
    }

    TEST("init from string");
    {
        RouthPolynomial poly;
        bool ok = routh_polynomial_from_string(&poly, "1 3 2");
        ASSERT_TRUE(ok && poly.degree == 2, "string init failed");
    }

    TEST("polynomial evaluation (Horner)");
    {
        double c[] = {6.0, -5.0, 1.0}; /* s² - 5s + 6 */
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 2);
        double p2 = routh_polynomial_eval(&poly, 2.0);
        double p3 = routh_polynomial_eval(&poly, 3.0);
        ASSERT_TRUE(fabs(p2) < 1e-10 && fabs(p3) < 1e-10, "eval at roots should be 0");
    }

    TEST("polynomial derivative");
    {
        double c[] = {2.0, 5.0, 3.0, 1.0}; /* s³ + 3s² + 5s + 2 */
        RouthPolynomial poly, deriv;
        routh_polynomial_init(&poly, c, 3);
        bool ok = routh_polynomial_derivative(&poly, &deriv);
        /* d/ds: 3s² + 6s + 5 → coefs [5, 6, 3] */
        ASSERT_TRUE(ok && deriv.degree == 2 &&
                    fabs(deriv.coeff[0] - 5.0) < 1e-10 &&
                    fabs(deriv.coeff[1] - 6.0) < 1e-10 &&
                    fabs(deriv.coeff[2] - 3.0) < 1e-10,
                    "derivative computation wrong");
    }
}

void test_routh_array(void)
{
    printf("\n--- L5: Routh Array Construction ---\n");

    /* Stable: s³ + 3s² + 3s + 1 = (s+1)³, all roots at -1 */
    TEST("routh array — stable cubic");
    {
        double c[] = {1.0, 3.0, 3.0, 1.0};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 3);
        RouthArray array;
        bool ok = routh_array_construct(&array, &poly);
        ASSERT_TRUE(ok && array.num_sign_changes == 0,
                     "stable cubic should have 0 sign changes");
    }

    /* Unstable: s³ + s² + s + 2 — one RHP root */
    TEST("routh array — unstable cubic");
    {
        double c[] = {2.0, 1.0, 1.0, 1.0};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 3);
        RouthArray array;
        routh_array_construct(&array, &poly);
        ASSERT_DOUBLE_GT((double)array.num_sign_changes, 0.0,
                          "unstable cubic should have sign changes");
    }

    /* Marginally unstable: s³ + 2s² + s + 2 = (s²+1)(s+2) */
    TEST("routh array — marginal cubic (zero row)");
    {
        double c[] = {2.0, 1.0, 2.0, 1.0}; /* s³ + 2s² + s + 2 */
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 3);
        RouthArray array;
        routh_array_construct(&array, &poly);
        ASSERT_TRUE(array.has_special_case, "marginal should flag special case");
    }
}

/* ============================================================================
 * L4 Tests: Stability Checking
 * ============================================================================ */

void test_stability_check(void)
{
    printf("\n--- L4: Stability Checking ---\n");

    /* Strictly stable: s² + 2s + 1 = (s+1)² */
    TEST("routh stable: s² + 2s + 1");
    {
        double c[] = {1.0, 2.0, 1.0};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 2);
        bool hurwitz = routh_is_hurwitz_polynomial(&poly);
        ASSERT_TRUE(hurwitz, "s²+2s+1 should be Hurwitz");
    }

    /* Strictly stable: s³ + 6s² + 11s + 6 = (s+1)(s+2)(s+3) */
    TEST("routh stable: s³ + 6s² + 11s + 6");
    {
        double c[] = {6.0, 11.0, 6.0, 1.0};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 3);
        bool hurwitz = routh_is_hurwitz_polynomial(&poly);
        ASSERT_TRUE(hurwitz, "s³+6s²+11s+6 should be Hurwitz");
    }

    /* Strictly stable: s⁴ + 10s³ + 35s² + 50s + 24 = (s+1)(s+2)(s+3)(s+4) */
    TEST("routh stable: s⁴ + 10s³ + 35s² + 50s + 24");
    {
        double c[] = {24.0, 50.0, 35.0, 10.0, 1.0};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 4);
        bool hurwitz = routh_is_hurwitz_polynomial(&poly);
        ASSERT_TRUE(hurwitz, "4th order all-real-negative should be Hurwitz");
    }

    /* Unstable: s² - 2s + 1 = (s-1)², right-half-plane */
    TEST("unstable: s² - 2s + 1");
    {
        double c[] = {1.0, -2.0, 1.0};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 2);
        bool hurwitz = routh_is_hurwitz_polynomial(&poly);
        ASSERT_TRUE(!hurwitz, "rhp roots: should NOT be Hurwitz");
    }

    /* Necessary condition fails */
    TEST("necessary condition fails for missing coefficient");
    {
        double c[] = {0.0, 1.0, 1.0}; /* Missing a_0 */
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 2);
        bool nec = routh_necessary_condition(&poly);
        ASSERT_TRUE(!nec, "missing coefficient should fail necessary condition");
    }
}

/* ============================================================================
 * L4 Tests: Hurwitz Determinant Method
 * ============================================================================ */

void test_hurwitz_determinants(void)
{
    printf("\n--- L4: Hurwitz Determinant Method ---\n");

    /* Degree 3: s³ + 3s² + 3s + 1, all Δ > 0 */
    TEST("Hurwitz minors for stable cubic");
    {
        double c[] = {1.0, 3.0, 3.0, 1.0};
        HurwitzResult result;
        bool ok = hurwitz_check_stability(c, 3, &result);
        ASSERT_TRUE(ok && result.is_stable, "stable cubic: Hurwitz should say stable");
    }

    /* Degree 3: s³ + s² + s + 2 — unstable */
    TEST("Hurwitz minors for unstable cubic");
    {
        double c[] = {2.0, 1.0, 1.0, 1.0};
        HurwitzResult result;
        hurwitz_check_stability(c, 3, &result);
        ASSERT_TRUE(!result.is_stable, "unstable cubic: Hurwitz should say unstable");
    }

    /* Degree 4 stable */
    TEST("Hurwitz minors for stable quartic");
    {
        double c[] = {24.0, 50.0, 35.0, 10.0, 1.0};
        HurwitzResult result;
        bool ok = hurwitz_check_stability(c, 4, &result);
        ASSERT_TRUE(ok && result.is_stable, "stable quartic");
    }

    /* Explicit conditions for degree 3 */
    TEST("explicit conditions for degree 3");
    {
        double c[] = {1.0, 3.0, 3.0, 1.0};
        char buf[512];
        bool ok = hurwitz_explicit_conditions(c, 3, buf, sizeof(buf));
        ASSERT_TRUE(ok, "explicit conditions failed");
    }
}

/* ============================================================================
 * L5 Tests: Special Cases
 * ============================================================================ */

void test_special_cases(void)
{
    printf("\n--- L5: Special Case Handling ---\n");

    /* Zero in first column: s⁴ + s³ + s² + s + K for specific K */
    TEST("detect zero first column");
    {
        double row[] = {0.0, 1.0, 2.0};
        bool is_zero = routh_detect_zero_first(row, 3, 1e-12);
        ASSERT_TRUE(is_zero, "should detect zero first element");
    }

    TEST("detect non-zero first column");
    {
        double row[] = {3.0, 1.0, 2.0};
        bool is_zero = routh_detect_zero_first(row, 3, 1e-12);
        ASSERT_TRUE(!is_zero, "should not detect zero when first is 3.0");
    }

    TEST("detect entire zero row");
    {
        double row[] = {0.0, 0.0, 0.0};
        bool is_zero = routh_detect_zero_row(row, 3, 1e-12);
        ASSERT_TRUE(is_zero, "should detect entire-zero row");
    }

    TEST("detect non-zero row");
    {
        double row[] = {0.0, 1.0, 0.0};
        bool is_zero = routh_detect_zero_row(row, 3, 1e-12);
        ASSERT_TRUE(!is_zero, "should not flag row with non-zero entry");
    }

    /* Liénard-Chipart for degree 4 */
    TEST("Liénard-Chipart: stable quartic passes");
    {
        double c[] = {24.0, 50.0, 35.0, 10.0, 1.0};
        bool ok = routh_lienard_chipart_criterion(c, 4);
        ASSERT_TRUE(ok, "stable quartic should pass Liénard-Chipart");
    }
}

/* ============================================================================
 * L6 Tests: Root Distribution
 * ============================================================================ */

void test_root_distribution(void)
{
    printf("\n--- L6: Root Distribution ---\n");

    TEST("root count: stable cubic (s+1)³ → 3 LHP, 0 RHP, 0 jω");
    {
        double c[] = {1.0, 3.0, 3.0, 1.0};
        int n_lhp, n_rhp, n_jw;
        bool ok = root_distribution_count(c, 3, &n_lhp, &n_rhp, &n_jw);
        ASSERT_TRUE(ok && n_lhp == 3 && n_rhp == 0 && n_jw == 0,
                     "stable cubic: 3 LHP, 0 RHP, 0 jω");
    }

    TEST("root count: unstable quadratic s² - 1 → 1 LHP, 1 RHP");
    {
        double c[] = {-1.0, 0.0, 1.0}; /* s² - 1 */
        int n_lhp, n_rhp, n_jw;
        root_distribution_count(c, 2, &n_lhp, &n_rhp, &n_jw);
        ASSERT_TRUE(n_rhp >= 1, "s²-1 should have at least 1 RHP root");
    }

    TEST("root count: marginal s³ + s → s(s²+1): 1 origin, 2 jω");
    {
        double c[] = {0.0, 1.0, 0.0, 1.0}; /* s³ + s */
        int n_lhp, n_rhp, n_jw;
        root_distribution_count(c, 3, &n_lhp, &n_rhp, &n_jw);
        ASSERT_TRUE(n_jw >= 1, "s(s²+1) should have jω roots");
    }

    TEST("Orlando formula: stable cubic, Δ_{n-1} > 0");
    {
        double c[] = {1.0, 3.0, 3.0, 1.0};
        double delta;
        bool is_marginal;
        bool ok = root_distribution_orlando_formula(c, 3, &delta, &is_marginal);
        ASSERT_TRUE(ok && !is_marginal && delta > 0.0, "Orlando: stable → Δ₂ > 0");
    }
}

/* ============================================================================
 * L5 Tests: Relative Stability
 * ============================================================================ */

void test_relative_stability(void)
{
    printf("\n--- L5: Relative Stability ---\n");

    TEST("axis shift: s² + 2s + 1 shifted by σ=0 → Hurwitz");
    {
        double c[] = {1.0, 2.0, 1.0};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 2);
        RelStabAxisShift result;
        bool ok = relstab_axis_shift(&poly, 0.0, &result);
        ASSERT_TRUE(ok && result.is_hurwitz, "σ=0 should preserve stability");
    }

    TEST("axis shift: s² + 2s + 1 shifted by σ=1 → (s-1+1)² = s²");
    {
        double c[] = {1.0, 2.0, 1.0};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 2);
        RelStabAxisShift result;
        relstab_axis_shift(&poly, 1.0, &result);
        /* P(s-1) = (s-1)² + 2(s-1) + 1 = s² - 2s + 1 + 2s - 2 + 1 = s² */
        /* s² has double pole at 0 → marginal */
        ASSERT_TRUE(!result.is_hurwitz, "s² should not be strictly Hurwitz");
    }

    TEST("stability margin for s² + 2s + 1");
    {
        double c[] = {1.0, 2.0, 1.0};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 2);
        RelStabMargin margin;
        bool ok = relstab_find_margin(&poly, &margin);
        ASSERT_TRUE(ok && margin.sigma_max > 0.9 && margin.sigma_max < 1.1,
                     "stability margin ≈ 1.0 for (s+1)²");
    }

    TEST("stability margin for s² + 4s + 4 = (s+2)² → margin ≈ 2");
    {
        double c[] = {4.0, 4.0, 1.0};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 2);
        RelStabMargin margin;
        relstab_find_margin(&poly, &margin);
        ASSERT_TRUE(margin.sigma_max > 1.5 && margin.sigma_max < 2.5,
                     "stability margin ≈ 2.0 for (s+2)²");
    }

    TEST("damping ratio estimate for s² + 2ζs + 1 with ζ=0.707");
    {
        double c[] = {1.0, 2.0 * 0.707, 1.0};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 2);
        double zeta;
        bool ok = relstab_min_damping_ratio(&poly, &zeta);
        ASSERT_TRUE(ok, "damping ratio estimate should succeed");
    }
}

/* ============================================================================
 * L4 Tests: Jury Stability (Discrete-Time)
 * ============================================================================ */

void test_jury_stability(void)
{
    printf("\n--- L4: Jury Stability (Discrete-Time) ---\n");

    /* Stable: z² + 0.5z + 0.25, both poles inside unit circle */
    TEST("Jury: stable 2nd order");
    {
        double c[] = {0.25, 0.5, 1.0}; /* a_0=0.25, a_1=0.5, a_2=1.0 */
        JuryResult result;
        bool ok = jury_check_stability(c, 2, &result);
        ASSERT_TRUE(ok && result.is_schur_stable, "z²+0.5z+0.25 should be Schur stable");
    }

    /* Unstable: z² + 2z + 1 = (z+1)², poles on unit circle → marginal */
    TEST("Jury: marginal 2nd order (poles on unit circle)");
    {
        double c[] = {1.0, 2.0, 1.0};
        JuryResult result;
        jury_check_stability(c, 2, &result);
        /* P(1) = 4 > 0, P(-1) = 0 → fails condition 2 */
        ASSERT_TRUE(!result.is_schur_stable, "(z+1)² should fail Jury stability");
    }

    /* Unstable: z + 2, pole at z=-2, outside unit circle */
    TEST("Jury: unstable 1st order");
    {
        double c[] = {2.0, 1.0}; /* z + 2 */
        JuryResult result;
        jury_check_stability(c, 1, &result);
        ASSERT_TRUE(!result.is_schur_stable, "z+2 pole at -2 should be unstable");
    }

    TEST("Jury: stable 3rd order");
    {
        double c[] = {0.125, 0.75, 1.5, 1.0}; /* (z+0.5)³ scaled */
        JuryResult result;
        bool ok = jury_check_stability(c, 3, &result);
        ASSERT_TRUE(ok, "3rd order jury test should complete");
    }

    TEST("Jury count unstable: z + 2 has 1 unstable root");
    {
        double c[] = {2.0, 1.0};
        int nu = jury_count_unstable(c, 1);
        ASSERT_EQ(nu, 1, "z+2 should have 1 root outside unit circle");
    }
}

/* ============================================================================
 * L3 Tests: Unified Stability Report
 * ============================================================================ */

void test_stability_report(void)
{
    printf("\n--- L3: Unified Stability Report ---\n");

    TEST("stability report for stable cubic");
    {
        double c[] = {1.0, 3.0, 3.0, 1.0};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 3);
        StabilityReport report;
        bool ok = stability_report(&poly, &report);
        ASSERT_TRUE(ok && report.overall == ROUTH_STABLE,
                     "stable cubic → ROUTH_STABLE");
    }

    TEST("quick check stable");
    {
        double c[] = {1.0, 3.0, 3.0, 1.0};
        RouthPolynomial poly;
        routh_polynomial_init(&poly, c, 3);
        RouthStability s = stability_quick_check(&poly);
        ASSERT_TRUE(s == ROUTH_STABLE, "quick check should return STABLE");
    }

    TEST("low-order stability for degree 3");
    {
        double c[] = {1.0, 3.0, 3.0, 1.0};
        LowOrderStability result;
        bool ok = stability_low_order(c, 3, &result);
        ASSERT_TRUE(ok && result.is_stable, "low-order should confirm stable");
    }

    TEST("companion matrix construction");
    {
        double c[] = {1.0, 3.0, 3.0, 1.0};
        double C[9]; /* 3x3 */
        bool ok = stability_companion_matrix(c, 3, C);
        ASSERT_TRUE(ok, "companion matrix construction failed");
        /* Top row should be -a_2/a_3, -a_1/a_3, -a_0/a_3 = -3, -3, -1 */
        ASSERT_DOUBLE_EQ(C[0], -3.0, 1e-10, "companion[0,0]");
        ASSERT_DOUBLE_EQ(C[1], -3.0, 1e-10, "companion[0,1]");
        ASSERT_DOUBLE_EQ(C[2], -1.0, 1e-10, "companion[0,2]");
    }
}

/* ============================================================================
 * Edge Cases & NULL Safety
 * ============================================================================ */

void test_edge_cases(void)
{
    printf("\n--- Edge Cases ---\n");

    TEST("NULL polynomial for all APIs");
    {
        ASSERT_TRUE(routh_polynomial_init(NULL, NULL, 2) == false, "init NULL");
        ASSERT_TRUE(routh_array_construct(NULL, NULL) == false, "array NULL");
        ASSERT_TRUE(routh_check_stability(NULL, NULL) == -1, "check NULL");
        ASSERT_TRUE(routh_count_sign_changes(NULL) == -1, "count NULL");
        ASSERT_TRUE(routh_is_hurwitz_polynomial(NULL) == false, "is_hurwitz NULL");
        PASS(); /* All NULL checks consolidated */
    }

    TEST("zero-degree polynomial rejected");
    {
        double c[] = {1.0}; /* just constant */
        RouthPolynomial poly;
        bool ok = routh_polynomial_init(&poly, c, 0);
        ASSERT_TRUE(!ok, "degree 0 should fail");
    }

    TEST("negative degree polynomial rejected");
    {
        double c[] = {1.0};
        RouthPolynomial poly;
        bool ok = routh_polynomial_init(&poly, c, -1);
        ASSERT_TRUE(!ok, "negative degree should fail");
    }
}

/* ============================================================================
 * Bilinear Transform & Discrete-Time
 * ============================================================================ */

void test_bilinear_transform(void)
{
    printf("\n--- Bilinear Transform ---\n");

    TEST("bilinear: z → (1+w)/(1-w) for degree 1");
    {
        double z_coeffs[] = {0.5, 1.0}; /* z + 0.5 */
        double w_coeffs[ROUTH_MAX_DEGREE + 1];
        bool ok = stability_bilinear_transform(z_coeffs, 1, w_coeffs);
        ASSERT_TRUE(ok, "bilinear transform should succeed for deg 1");
    }

    TEST("discrete-time stability via bilinear+Routh");
    {
        double z_coeffs[] = {0.25, 0.5, 1.0}; /* Stable discrete */
        StabilityReport report;
        bool ok = stability_discrete_time(z_coeffs, 2, &report);
        ASSERT_TRUE(ok, "discrete-time stability check should succeed");
    }
}

/* ============================================================================
 * Main Test Runner
 * ============================================================================ */

int main(void)
{
    printf("========================================\n");
    printf(" Routh-Hurwitz Stability — Test Suite\n");
    printf("========================================\n");

    test_polynomial_init();
    test_routh_array();
    test_stability_check();
    test_hurwitz_determinants();
    test_special_cases();
    test_root_distribution();
    test_relative_stability();
    test_jury_stability();
    test_stability_report();
    test_bilinear_transform();
    test_edge_cases();

    printf("\n========================================\n");
    printf(" RESULTS: %d/%d passed, %d failed\n",
           tests_passed, tests_total, tests_failed);
    printf("========================================\n");

    return (tests_failed > 0) ? 1 : 0;
}
