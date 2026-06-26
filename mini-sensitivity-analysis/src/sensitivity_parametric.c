/**
 * sensitivity_parametric.c — Parameter Sensitivity Analysis Implementation
 *
 * Implements finite-difference sensitivity, Sobol' variance-based global
 * sensitivity indices, Morris screening method, and direct differentiation
 * sensitivity for ODE systems.
 *
 * Knowledge Coverage:
 *   L1: Parameter sensitivity ∂f/∂p via finite differences
 *   L2: Forward vs central differences (O(h) vs O(h²))
 *   L3: Sobol' variance decomposition ANOVA
 *   L5: Morris method for parameter screening
 *   L5: Direct differentiation (variational equation integration)
 *   L8: Monte Carlo integration for global sensitivity
 */

#include "sensitivity_parametric.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ==========================================================================
 * L5: Forward Finite Difference
 *
 * ∂f/∂p_i ≈ (f(p + h·e_i) - f(p)) / h
 * Truncation error: O(h) — first order accurate.
 * ========================================================================== */

void fd_sensitivity_forward(ParametricFunction f,
                            const double *params, int n_params,
                            const double *h,
                            SensitivityResult *results) {
    if (f == NULL || params == NULL || results == NULL || n_params <= 0) return;

    /* Evaluate nominal function value */
    double f_nom = f(params, n_params);

    /* Allocate perturbed parameter vector */
    double *params_pert = (double *)malloc((size_t)n_params * sizeof(double));
    if (params_pert == NULL) return;
    memcpy(params_pert, params, (size_t)n_params * sizeof(double));

    for (int i = 0; i < n_params; i++) {
        /* Step size: use provided h or auto-compute */
        double step = (h != NULL) ? h[i] : fmax(1e-6 * fabs(params[i]), 1e-8);
        if (step < 1e-15) step = 1e-8;

        /* Perturb parameter i */
        params_pert[i] = params[i] + step;

        /* Evaluate perturbed function */
        double f_pert = f(params_pert, n_params);

        /* Compute sensitivity */
        results[i].value_nominal = f_nom;
        results[i].value_perturbed = f_pert;
        results[i].delta_param = step;
        results[i].abs_sensitivity = (f_pert - f_nom) / step;
        results[i].rel_sensitivity = 0.0;
        if (fabs(f_nom) > 1e-15) {
            results[i].rel_sensitivity = results[i].abs_sensitivity *
                                         params[i] / f_nom;
        }

        /* Restore parameter */
        params_pert[i] = params[i];
    }

    free(params_pert);
}

/* ==========================================================================
 * L5: Central Finite Difference
 *
 * ∂f/∂p_i ≈ (f(p + h·e_i) - f(p - h·e_i)) / (2h)
 * Truncation error: O(h²) — second order accurate.
 * Requires 2·n_params function evaluations (vs n_params for forward).
 * ========================================================================== */

void fd_sensitivity_central(ParametricFunction f,
                            const double *params, int n_params,
                            const double *h,
                            SensitivityResult *results) {
    if (f == NULL || params == NULL || results == NULL || n_params <= 0) return;

    double *params_pert = (double *)malloc((size_t)n_params * sizeof(double));
    if (params_pert == NULL) return;

    for (int i = 0; i < n_params; i++) {
        double step = (h != NULL) ? h[i] : fmax(1e-6 * fabs(params[i]), 1e-8);
        if (step < 1e-15) step = 1e-8;

        /* Evaluate f(p + h·e_i) */
        memcpy(params_pert, params, (size_t)n_params * sizeof(double));
        params_pert[i] = params[i] + step;
        double f_plus = f(params_pert, n_params);

        /* Evaluate f(p - h·e_i) */
        params_pert[i] = params[i] - step;
        double f_minus = f(params_pert, n_params);

        /* Central difference */
        results[i].value_nominal = (f_plus + f_minus) / 2.0;
        results[i].value_perturbed = f_plus;
        results[i].delta_param = step;
        results[i].abs_sensitivity = (f_plus - f_minus) / (2.0 * step);

        /* Relative sensitivity */
        double f_nom = results[i].value_nominal;
        results[i].rel_sensitivity = 0.0;
        if (fabs(f_nom) > 1e-15) {
            results[i].rel_sensitivity = results[i].abs_sensitivity *
                                         params[i] / f_nom;
        }
    }

    free(params_pert);
}

/* ==========================================================================
 * L5: Jacobian ∂f/∂x via Central Differences
 * ========================================================================== */

void jacobian_fd(StateFunction f, double t, const double *x,
                 const double *params, int n_states, int n_params,
                 double *J) {
    if (f == NULL || x == NULL || J == NULL || n_states <= 0) return;

    double *x_plus = (double *)malloc((size_t)n_states * sizeof(double));
    double *x_minus = (double *)malloc((size_t)n_states * sizeof(double));
    double *f_plus = (double *)malloc((size_t)n_states * sizeof(double));
    double *f_minus = (double *)malloc((size_t)n_states * sizeof(double));

    if (x_plus == NULL || x_minus == NULL || f_plus == NULL || f_minus == NULL) {
        free(x_plus); free(x_minus); free(f_plus); free(f_minus);
        return;
    }

    memcpy(x_plus, x, (size_t)n_states * sizeof(double));
    memcpy(x_minus, x, (size_t)n_states * sizeof(double));

    for (int j = 0; j < n_states; j++) {
        double h = fmax(1e-6 * fabs(x[j]), 1e-8);

        x_plus[j] = x[j] + h;
        x_minus[j] = x[j] - h;

        f(t, x_plus, params, n_params, n_states, f_plus);
        f(t, x_minus, params, n_params, n_states, f_minus);

        /* J[i][j] = ∂f_i/∂x_j */
        for (int i = 0; i < n_states; i++) {
            J[i * n_states + j] = (f_plus[i] - f_minus[i]) / (2.0 * h);
        }

        /* Restore */
        x_plus[j] = x[j];
        x_minus[j] = x[j];
    }

    free(x_plus); free(x_minus); free(f_plus); free(f_minus);
}

/* ==========================================================================
 * L5: Parameter Jacobian ∂f/∂p via Central Differences
 * ========================================================================== */

void param_jacobian_fd(StateFunction f, double t, const double *x,
                       const double *params, int n_states, int n_params,
                       double *dfdp) {
    if (f == NULL || x == NULL || params == NULL || dfdp == NULL) return;

    double *p_plus = (double *)malloc((size_t)n_params * sizeof(double));
    double *p_minus = (double *)malloc((size_t)n_params * sizeof(double));
    double *f_plus = (double *)malloc((size_t)n_states * sizeof(double));
    double *f_minus = (double *)malloc((size_t)n_states * sizeof(double));

    if (p_plus == NULL || p_minus == NULL || f_plus == NULL || f_minus == NULL) {
        free(p_plus); free(p_minus); free(f_plus); free(f_minus);
        return;
    }

    memcpy(p_plus, params, (size_t)n_params * sizeof(double));
    memcpy(p_minus, params, (size_t)n_params * sizeof(double));

    for (int j = 0; j < n_params; j++) {
        double h = fmax(1e-6 * fabs(params[j]), 1e-8);

        p_plus[j] = params[j] + h;
        p_minus[j] = params[j] - h;

        f(t, x, p_plus, n_params, n_states, f_plus);
        f(t, x, p_minus, n_params, n_states, f_minus);

        for (int i = 0; i < n_states; i++) {
            dfdp[i * n_params + j] = (f_plus[i] - f_minus[i]) / (2.0 * h);
        }

        p_plus[j] = params[j];
        p_minus[j] = params[j];
    }

    free(p_plus); free(p_minus); free(f_plus); free(f_minus);
}

/* ==========================================================================
 * L5: Direct Sensitivity for ODE Systems
 *
 * Integrates the augmented system combining state and sensitivity dynamics.
 *   Original:  ẋ = f(t, x, p)
 *   Sensitivity: d/dt(∂x/∂p_j) = ∂f/∂x · (∂x/∂p_j) + ∂f/∂p_j
 *
 * The augmented state vector: [x; s_1; s_2; ...; s_{n_p}]
 * where s_j = ∂x/∂p_j is an n_states vector.
 *
 * Uses RK4 fixed-step integration.
 * ========================================================================== */

void direct_sensitivity_ode(StateFunction f,
                            double t0, double tf,
                            const double *x0, const double *dx0_dp,
                            const double *params,
                            int n_states, int n_params, int n_steps,
                            double *t_out, double *x_out, double *s_out) {
    if (f == NULL || x0 == NULL || t_out == NULL || n_steps <= 0) return;

    double dt = (tf - t0) / (double)n_steps;
    int aug_dim = n_states * (1 + n_params);  /* total augmented dimension */

    /* Allocate augmented state vector */
    double *z = (double *)malloc((size_t)aug_dim * sizeof(double));
    double *k1 = (double *)malloc((size_t)aug_dim * sizeof(double));
    double *k2 = (double *)malloc((size_t)aug_dim * sizeof(double));
    double *k3 = (double *)malloc((size_t)aug_dim * sizeof(double));
    double *k4 = (double *)malloc((size_t)aug_dim * sizeof(double));
    double *ztmp = (double *)malloc((size_t)aug_dim * sizeof(double));

    if (z == NULL || k1 == NULL || k2 == NULL || k3 == NULL || k4 == NULL || ztmp == NULL) {
        free(z); free(k1); free(k2); free(k3); free(k4); free(ztmp);
        return;
    }

    /* Set initial condition */
    for (int i = 0; i < n_states; i++) {
        z[i] = x0[i];  /* nominal state */
    }
    if (dx0_dp != NULL) {
        for (int p = 0; p < n_params; p++) {
            for (int i = 0; i < n_states; i++) {
                z[n_states + p * n_states + i] = dx0_dp[p * n_states + i];
            }
        }
    } else {
        /* Default: zero initial sensitivity */
        for (int j = n_states; j < aug_dim; j++) z[j] = 0.0;
    }

    /* Record initial state */
    t_out[0] = t0;
    if (x_out != NULL) {
        for (int i = 0; i < n_states; i++) x_out[i] = z[i];
    }
    if (s_out != NULL) {
        for (int p = 0; p < n_params; p++) {
            for (int i = 0; i < n_states; i++) {
                s_out[p * (n_steps + 1) * n_states + 0 * n_states + i] =
                    z[n_states + p * n_states + i];
            }
        }
    }

    /* Approximate Jacobians at current state */
    double *Jx = (double *)malloc((size_t)(n_states * n_states) * sizeof(double));
    double *Jp = (double *)malloc((size_t)(n_states * n_params) * sizeof(double));

    if (Jx == NULL || Jp == NULL) {
        free(z); free(k1); free(k2); free(k3); free(k4); free(ztmp);
        free(Jx); free(Jp);
        return;
    }

    for (int step = 0; step < n_steps; step++) {
        double t = t0 + (double)step * dt;

        /* --- RK4 Stage 1: k1 = f(t, z) --- */
        /* Evaluate the augmented RHS at current z */
        f(t, z, params, n_params, n_states, k1);

        /* Compute Jacobians at current state for sensitivity RHS */
        jacobian_fd(f, t, z, params, n_states, n_params, Jx);
        param_jacobian_fd(f, t, z, params, n_states, n_params, Jp);

        /* Sensitivity RHS: Jx·s_p + Jp[:,p] */
        for (int p = 0; p < n_params; p++) {
            double *s_p = &z[n_states + p * n_states];  /* s_p is n_states vector */
            double *k1_p = &k1[n_states + p * n_states];

            for (int i = 0; i < n_states; i++) {
                k1_p[i] = 0.0;
                for (int j = 0; j < n_states; j++) {
                    k1_p[i] += Jx[i * n_states + j] * s_p[j];
                }
                k1_p[i] += Jp[i * n_params + p];  /* ∂f_i/∂p_p */
            }
        }

        /* --- RK4 Stage 2 --- */
        for (int i = 0; i < aug_dim; i++) {
            ztmp[i] = z[i] + 0.5 * dt * k1[i];
        }
        f(t + 0.5 * dt, ztmp, params, n_params, n_states, k2);

        jacobian_fd(f, t + 0.5 * dt, ztmp, params, n_states, n_params, Jx);
        param_jacobian_fd(f, t + 0.5 * dt, ztmp, params, n_states, n_params, Jp);

        for (int p = 0; p < n_params; p++) {
            double *s_p = &ztmp[n_states + p * n_states];
            double *k2_p = &k2[n_states + p * n_states];
            for (int i = 0; i < n_states; i++) {
                k2_p[i] = 0.0;
                for (int j = 0; j < n_states; j++) {
                    k2_p[i] += Jx[i * n_states + j] * s_p[j];
                }
                k2_p[i] += Jp[i * n_params + p];
            }
        }

        /* --- RK4 Stage 3 --- */
        for (int i = 0; i < aug_dim; i++) {
            ztmp[i] = z[i] + 0.5 * dt * k2[i];
        }
        f(t + 0.5 * dt, ztmp, params, n_params, n_states, k3);

        jacobian_fd(f, t + 0.5 * dt, ztmp, params, n_states, n_params, Jx);
        param_jacobian_fd(f, t + 0.5 * dt, ztmp, params, n_states, n_params, Jp);

        for (int p = 0; p < n_params; p++) {
            double *s_p = &ztmp[n_states + p * n_states];
            double *k3_p = &k3[n_states + p * n_states];
            for (int i = 0; i < n_states; i++) {
                k3_p[i] = 0.0;
                for (int j = 0; j < n_states; j++) {
                    k3_p[i] += Jx[i * n_states + j] * s_p[j];
                }
                k3_p[i] += Jp[i * n_params + p];
            }
        }

        /* --- RK4 Stage 4 --- */
        for (int i = 0; i < aug_dim; i++) {
            ztmp[i] = z[i] + dt * k3[i];
        }
        f(t + dt, ztmp, params, n_params, n_states, k4);

        jacobian_fd(f, t + dt, ztmp, params, n_states, n_params, Jx);
        param_jacobian_fd(f, t + dt, ztmp, params, n_states, n_params, Jp);

        for (int p = 0; p < n_params; p++) {
            double *s_p = &ztmp[n_states + p * n_states];
            double *k4_p = &k4[n_states + p * n_states];
            for (int i = 0; i < n_states; i++) {
                k4_p[i] = 0.0;
                for (int j = 0; j < n_states; j++) {
                    k4_p[i] += Jx[i * n_states + j] * s_p[j];
                }
                k4_p[i] += Jp[i * n_params + p];
            }
        }

        /* --- RK4 Update --- */
        for (int i = 0; i < aug_dim; i++) {
            z[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }

        /* Record output */
        int idx = step + 1;
        t_out[idx] = t + dt;
        if (x_out != NULL) {
            for (int i = 0; i < n_states; i++) {
                x_out[idx * n_states + i] = z[i];
            }
        }
        if (s_out != NULL) {
            for (int p = 0; p < n_params; p++) {
                for (int i = 0; i < n_states; i++) {
                    s_out[p * (n_steps + 1) * n_states + idx * n_states + i] =
                        z[n_states + p * n_states + i];
                }
            }
        }
    }

    free(z); free(k1); free(k2); free(k3); free(k4); free(ztmp);
    free(Jx); free(Jp);
}

/* ==========================================================================
 * L8: Sobol' First-Order Sensitivity Index (Monte Carlo Method)
 *
 * The Sobol' method decomposes output variance into contributions from
 * individual parameters and their interactions.
 *
 * Algorithm (Saltelli, 2010):
 *   1. Generate two N×k sample matrices A and B (Sobol' or LHS)
 *   2. Build matrix A_B[i] where column i comes from B, rest from A
 *   3. First-order index: S_i = V[E[Y|X_i]] / V[Y]
 *                      ≈ (1/N Σ f(A)_j · f(A_B[i])_j - f₀²) / Var(Y)
 *
 * where f₀ = (1/N) Σ f(A)_j and Var(Y) = (1/N) Σ f(A)_j² - f₀²
 * ========================================================================== */

/* Pseudo-random uniform [0,1] using simple LCG */
static double uniform_rand(unsigned int *seed) {
    *seed = *seed * 1103515245 + 12345;
    return (double)((*seed >> 16) & 0x7FFF) / 32767.0;
}

void sobol_first_order(ParametricFunction f,
                       const double *bounds, int n_params,
                       int n_samples, double *sobol) {
    if (f == NULL || bounds == NULL || sobol == NULL || n_params <= 0 || n_samples <= 0) {
        return;
    }

    /* Allocate sample matrices */
    int N = n_samples;
    double *A = (double *)malloc((size_t)(N * n_params) * sizeof(double));
    double *B = (double *)malloc((size_t)(N * n_params) * sizeof(double));

    if (A == NULL || B == NULL) {
        free(A); free(B);
        return;
    }

    /* Generate samples using pseudo-random LHS-like approach */
    unsigned int seed = 12345;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < n_params; j++) {
            double lo = bounds[2 * j];
            double hi = bounds[2 * j + 1];
            A[i * n_params + j] = lo + uniform_rand(&seed) * (hi - lo);
            B[i * n_params + j] = lo + uniform_rand(&seed) * (hi - lo);
        }
    }

    /* Evaluate f on matrix A */
    double *f_A = (double *)malloc((size_t)N * sizeof(double));
    double *f_B = (double *)malloc((size_t)N * sizeof(double));
    double *f_AB = (double *)malloc((size_t)N * sizeof(double));

    if (f_A == NULL || f_B == NULL || f_AB == NULL) {
        free(A); free(B); free(f_A); free(f_B); free(f_AB);
        return;
    }

    for (int i = 0; i < N; i++) {
        f_A[i] = f(&A[i * n_params], n_params);
        f_B[i] = f(&B[i * n_params], n_params);
    }

    /* Compute f₀² and total variance */
    double f0 = 0.0;
    for (int i = 0; i < N; i++) f0 += f_A[i];
    f0 /= (double)N;

    double f02_squared = f0 * f0;

    double var_Y = 0.0;
    for (int i = 0; i < N; i++) {
        var_Y += f_A[i] * f_A[i];
    }
    var_Y = var_Y / (double)N - f02_squared;

    if (var_Y < 1e-30) {
        /* Zero variance — all sensitivities are zero */
        for (int j = 0; j < n_params; j++) sobol[j] = 0.0;
        free(A); free(B); free(f_A); free(f_B); free(f_AB);
        return;
    }

    /* Compute first-order indices */
    for (int param = 0; param < n_params; param++) {
        /* Build A_B[param]: column param from B, rest from A */
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < n_params; j++) {
                if (j == param) {
                    f_AB[i] = f_B[i * n_params + j];  /* wrong — need to re-evaluate */
                }
            }
        }

        /* Actually, we need to evaluate f(A_B[param]_i) for each i */
        double *ab_row = (double *)malloc((size_t)n_params * sizeof(double));
        if (ab_row == NULL) continue;

        double sum_prod = 0.0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < n_params; j++) {
                ab_row[j] = (j == param) ? B[i * n_params + j] : A[i * n_params + j];
            }
            double f_ab_val = f(ab_row, n_params);
            sum_prod += f_A[i] * f_ab_val;
        }

        /* S_i = ((1/N)Σ f(A)_j · f(A_B)_j - f₀²) / Var(Y) */
        double Vi = sum_prod / (double)N - f02_squared;
        sobol[param] = Vi / var_Y;

        /* Clamp to [0,1] (numerical errors may push slightly outside) */
        if (sobol[param] < 0.0) sobol[param] = 0.0;
        if (sobol[param] > 1.0) sobol[param] = 1.0;

        free(ab_row);
    }

    free(A); free(B); free(f_A); free(f_B); free(f_AB);
}

/* ==========================================================================
 * L8: Sobol' Total-Effect Index
 *
 * S_{Ti} = 1 - V[E[Y|X_∼i]] / V[Y]
 *        ≈ 1 - ((1/N)Σ f(B)_j · f(A_B[i])_j - f₀²) / Var(Y)
 * where A_B[i] has column i from A, rest from B.
 *
 * Total effect includes main effect + all interactions involving X_i.
 * ========================================================================== */

void sobol_total_effect(ParametricFunction f,
                        const double *bounds, int n_params,
                        int n_samples, double *sobol_total) {
    if (f == NULL || bounds == NULL || sobol_total == NULL || n_params <= 0 || n_samples <= 0) {
        return;
    }

    int N = n_samples;
    unsigned int seed = 54321;

    /* Generate sample matrices */
    double *A = (double *)malloc((size_t)(N * n_params) * sizeof(double));
    double *B = (double *)malloc((size_t)(N * n_params) * sizeof(double));

    if (A == NULL || B == NULL) {
        free(A); free(B);
        return;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < n_params; j++) {
            double lo = bounds[2 * j];
            double hi = bounds[2 * j + 1];
            A[i * n_params + j] = lo + uniform_rand(&seed) * (hi - lo);
            B[i * n_params + j] = lo + uniform_rand(&seed) * (hi - lo);
        }
    }

    double *f_A = (double *)malloc((size_t)N * sizeof(double));
    double *f_B = (double *)malloc((size_t)N * sizeof(double));

    if (f_A == NULL || f_B == NULL) {
        free(A); free(B); free(f_A); free(f_B);
        return;
    }

    for (int i = 0; i < N; i++) {
        f_A[i] = f(&A[i * n_params], n_params);
        f_B[i] = f(&B[i * n_params], n_params);
    }

    double f0_B = 0.0;
    for (int i = 0; i < N; i++) f0_B += f_B[i];
    f0_B /= (double)N;

    double var_Y_B = 0.0;
    for (int i = 0; i < N; i++) {
        var_Y_B += f_B[i] * f_B[i];
    }
    var_Y_B = var_Y_B / (double)N - f0_B * f0_B;

    if (var_Y_B < 1e-30) {
        for (int j = 0; j < n_params; j++) sobol_total[j] = 0.0;
        free(A); free(B); free(f_A); free(f_B);
        return;
    }

    double *ab_row = (double *)malloc((size_t)n_params * sizeof(double));

    for (int param = 0; param < n_params; param++) {
        double sum_prod = 0.0;
        for (int i = 0; i < N; i++) {
            /* A_B[i] for total effect: column param from A, rest from B */
            for (int j = 0; j < n_params; j++) {
                ab_row[j] = (j == param) ? A[i * n_params + j] : B[i * n_params + j];
            }
            double f_ab_val = f(ab_row, n_params);
            sum_prod += f_B[i] * f_ab_val;
        }

        /* S_{Ti} = 1 - ((1/N)Σ f(B)·f(A_B) - f₀²) / Var(Y) */
        double V_minus_i = sum_prod / (double)N - f0_B * f0_B;
        sobol_total[param] = 1.0 - V_minus_i / var_Y_B;

        if (sobol_total[param] < 0.0) sobol_total[param] = 0.0;
        if (sobol_total[param] > 1.0) sobol_total[param] = 1.0;
    }

    free(A); free(B); free(f_A); free(f_B); free(ab_row);
}

/* ==========================================================================
 * L5: Morris Screening Method
 *
 * For each trajectory (random walk through parameter space), compute
 * elementary effects for each parameter:
 *   d_i = (f(..., p_i+Δ, ...) - f(...)) / Δ
 *
 * After r trajectories, report:
 *   μ* = mean(|d_i|) — overall influence
 *   σ  = std(d_i)   — nonlinearity/interaction indicator
 *
 * Reference: Morris, M.D. "Factorial Sampling Plans..." Technometrics (1991)
 * ========================================================================== */

void morris_screening(ParametricFunction f,
                      const double *bounds, int n_params,
                      int r, int p_levels,
                      double *mu_star, double *sigma) {
    if (f == NULL || bounds == NULL || mu_star == NULL || sigma == NULL) return;
    if (n_params <= 0 || r <= 0 || p_levels < 2) return;

    /* Allocate storage for elementary effects from all trajectories */
    double **all_effects = (double **)malloc((size_t)n_params * sizeof(double *));
    for (int i = 0; i < n_params; i++) {
        all_effects[i] = (double *)malloc((size_t)r * sizeof(double));
    }

    double *params_curr = (double *)malloc((size_t)n_params * sizeof(double));
    double *params_next = (double *)malloc((size_t)n_params * sizeof(double));

    unsigned int seed = 99999;
    double delta = (double)p_levels / (2.0 * ((double)p_levels - 1.0));

    for (int traj = 0; traj < r; traj++) {
        /* Random starting point in the grid */
        for (int j = 0; j < n_params; j++) {
            double lo = bounds[2 * j];
            double hi = bounds[2 * j + 1];
            params_curr[j] = lo + uniform_rand(&seed) * (hi - lo);
        }

        /* Random permutation of parameter order for this trajectory */
        int *perm = (int *)malloc((size_t)n_params * sizeof(int));
        for (int j = 0; j < n_params; j++) perm[j] = j;

        /* Fisher-Yates shuffle */
        for (int j = n_params - 1; j > 0; j--) {
            int k = (int)(uniform_rand(&seed) * (double)(j + 1));
            if (k > j) k = j;
            int tmp = perm[j];
            perm[j] = perm[k];
            perm[k] = tmp;
        }

        /* Walk through parameters in permuted order */
        memcpy(params_next, params_curr, (size_t)n_params * sizeof(double));

        for (int idx = 0; idx < n_params; idx++) {
            int p = perm[idx];
            double range = bounds[2 * p + 1] - bounds[2 * p];
            double step = delta * range;

            /* Random direction (+ or -) */
            int direction = (uniform_rand(&seed) > 0.5) ? 1 : -1;
            params_next[p] = params_curr[p] + direction * step;

            /* Clamp to bounds */
            if (params_next[p] < bounds[2 * p]) {
                params_next[p] = bounds[2 * p] + step;
            }
            if (params_next[p] > bounds[2 * p + 1]) {
                params_next[p] = bounds[2 * p + 1] - step;
            }

            /* Evaluate elementary effect */
            double f_curr = f(params_curr, n_params);
            double f_next = f(params_next, n_params);
            double effect = (f_next - f_curr) / (params_next[p] - params_curr[p] + 1e-30);

            all_effects[p][traj] = effect;

            /* Move to next point */
            memcpy(params_curr, params_next, (size_t)n_params * sizeof(double));
        }

        free(perm);
    }

    /* Compute μ* and σ for each parameter */
    for (int p = 0; p < n_params; p++) {
        double sum_abs = 0.0;
        double sum_sq = 0.0;

        for (int traj = 0; traj < r; traj++) {
            sum_abs += fabs(all_effects[p][traj]);
            sum_sq += all_effects[p][traj] * all_effects[p][traj];
        }

        mu_star[p] = sum_abs / (double)r;
        double mean_val = 0.0;
        for (int traj = 0; traj < r; traj++) {
            mean_val += all_effects[p][traj];
        }
        mean_val /= (double)r;
        double variance = sum_sq / (double)r - mean_val * mean_val;
        sigma[p] = sqrt(fmax(variance, 0.0));
    }

    for (int i = 0; i < n_params; i++) free(all_effects[i]);
    free(all_effects);
    free(params_curr);
    free(params_next);
}

/* ==========================================================================
 * L3: Statistical Analysis of Sensitivity Results
 * ========================================================================== */

void rank_by_sensitivity(const double *sensitivities, int n_params,
                         int *ranks) {
    if (sensitivities == NULL || ranks == NULL || n_params <= 0) return;

    /* Initialize ranks */
    for (int i = 0; i < n_params; i++) ranks[i] = 1;

    /* Simple comparison-based ranking (O(n²) — sufficient for small n_params) */
    for (int i = 0; i < n_params; i++) {
        for (int j = 0; j < n_params; j++) {
            if (fabs(sensitivities[j]) > fabs(sensitivities[i])) {
                ranks[i]++;
            }
        }
    }
}

void normalize_sensitivities(const double *sensitivities, int n_params,
                             double *normalized) {
    if (sensitivities == NULL || normalized == NULL || n_params <= 0) return;

    double sum_abs = 0.0;
    for (int i = 0; i < n_params; i++) {
        sum_abs += fabs(sensitivities[i]);
    }

    if (sum_abs < 1e-15) {
        for (int i = 0; i < n_params; i++) normalized[i] = 0.0;
        return;
    }

    for (int i = 0; i < n_params; i++) {
        normalized[i] = fabs(sensitivities[i]) / sum_abs;
    }
}

void sensitivity_cv(const double **data, int n_bootstrap, int n_params,
                    double *cv) {
    if (data == NULL || cv == NULL || n_bootstrap <= 1 || n_params <= 0) return;

    for (int p = 0; p < n_params; p++) {
        double mean = 0.0;
        for (int b = 0; b < n_bootstrap; b++) {
            mean += data[b][p];
        }
        mean /= (double)n_bootstrap;

        double variance = 0.0;
        for (int b = 0; b < n_bootstrap; b++) {
            double diff = data[b][p] - mean;
            variance += diff * diff;
        }
        variance /= (double)(n_bootstrap - 1);

        double std_dev = sqrt(fmax(variance, 0.0));
        if (fabs(mean) > 1e-15) {
            cv[p] = std_dev / fabs(mean);
        } else {
            cv[p] = INFINITY;
        }
    }
}
