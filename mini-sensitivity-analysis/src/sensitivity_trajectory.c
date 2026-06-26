/**
 * sensitivity_trajectory.c — Trajectory Sensitivity Analysis
 *
 * Implements trajectory sensitivity for ODE systems, including
 * periodic orbit sensitivity, initial condition sensitivity,
 * and discrete-time trajectory sensitivity.
 *
 * Knowledge Coverage:
 *   L2: Trajectory sensitivity equations
 *   L4: Variational equation for periodic orbits
 *   L5: Discrete-time trajectory sensitivity
 *   L5: Initial condition sensitivity via state transition matrix
 *   L8: Periodic orbit sensitivity via monodromy matrix
 */

#include "sensitivity_core.h"
#include "sensitivity_eigenvalue.h"
#include "sensitivity_time_domain.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* ==========================================================================
 * L4: State Transition Matrix Sensitivity
 *
 * The state transition matrix Φ(t, t₀) = ∂x(t)/∂x(t₀) satisfies:
 *   dΦ/dt = A(t)·Φ,  Φ(t₀, t₀) = I
 * where A(t) = ∂f/∂x evaluated along the nominal trajectory.
 *
 * The sensitivity of x(t) to initial condition x(0) = x₀ is:
 *   ∂x(t)/∂x₀ = Φ(t, 0)
 *
 * The sensitivity of x(t) to a parameter p is:
 *   ∂x(t)/∂p = ∫₀^t Φ(t, τ)·(∂f/∂p)(τ) dτ + Φ(t, 0)·(∂x₀/∂p)
 * ========================================================================== */

/**
 * Compute the state transition matrix Φ(t_f, 0) by integrating
 * the variational equation alongside the nominal trajectory.
 *
 * The augmented system has dimension n + n² (state + vectorized Φ).
 *
 * @param A_func function returning A(t) = ∂f/∂x at time t
 * @param n state dimension
 * @param tf final time
 * @param n_steps integration steps
 * @param Phi output state transition matrix (n×n, row-major)
 */
void state_transition_matrix(
    void (*A_func)(double t, const double *x, int n, double *A_mat),
    void (*dynamics)(double t, const double *x, int n, double *dxdt),
    int n, const double *x0, double tf, int n_steps,
    double *Phi) {

    if (A_func == NULL || dynamics == NULL || Phi == NULL || n <= 0) return;

    double dt = tf / (double)n_steps;

    /* Allocate state and Phi */
    double *x = (double *)malloc((size_t)n * sizeof(double));
    double *Phi_flat = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *A_mat = (double *)malloc((size_t)(n * n) * sizeof(double));

    if (x == NULL || Phi_flat == NULL || A_mat == NULL) {
        free(x); free(Phi_flat); free(A_mat);
        return;
    }

    /* Initialize */
    memcpy(x, x0, (size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Phi_flat[i * n + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    /* Forward Euler integration for both x and Φ */
    for (int step = 0; step < n_steps; step++) {
        double t = (double)step * dt;

        /* Compute A(t) */
        A_func(t, x, n, A_mat);

        /* Compute dx/dt */
        double *dxdt = (double *)malloc((size_t)n * sizeof(double));
        if (dxdt == NULL) break;
        dynamics(t, x, n, dxdt);

        /* Compute dΦ/dt = A(t)·Φ */
        double *dPhi = (double *)malloc((size_t)(n * n) * sizeof(double));
        if (dPhi == NULL) {
            free(dxdt);
            break;
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                dPhi[i * n + j] = 0.0;
                for (int k = 0; k < n; k++) {
                    dPhi[i * n + j] += A_mat[i * n + k] * Phi_flat[k * n + j];
                }
            }
        }

        /* Euler step */
        for (int i = 0; i < n; i++) {
            x[i] += dt * dxdt[i];
        }
        for (int i = 0; i < n * n; i++) {
            Phi_flat[i] += dt * dPhi[i];
        }

        free(dxdt);
        free(dPhi);
    }

    memcpy(Phi, Phi_flat, (size_t)(n * n) * sizeof(double));

    free(x); free(Phi_flat); free(A_mat);
}

/* ==========================================================================
 * L4: Initial Condition Sensitivity
 *
 * For ẋ = f(x, p), the sensitivity to initial condition x(0) = x₀:
 *   s_{x₀}(t) = ∂x(t)/∂x₀ = Φ(t, 0)
 *
 * This is the state transition matrix. For LTI systems:
 *   Φ(t, 0) = e^{A·t}
 * ========================================================================== */

/**
 * Compute the sensitivity of x(t) to initial condition perturbation.
 * Returns s_{ij}(t) = ∂x_i(t)/∂x₀_j for all i,j.
 *
 * For LTI systems, Φ(t) = e^{A·t}.
 *
 * @param A state matrix (n×n, row-major)
 * @param n state dimension
 * @param t time
 * @param s_IC output sensitivity matrix (n×n, row-major): s_IC[i*n+j] = ∂x_i/∂x₀_j
 */
void initial_condition_sensitivity_lti(const double *A, int n,
                                        double t, double *s_IC) {
    if (A == NULL || s_IC == NULL || n <= 0) return;

    /* s_IC = e^{A·t} = Φ(t) */
    matrix_exponential(A, n, t, s_IC);
}

/**
 * Trace-based measure of initial condition sensitivity:
 *   μ(t) = trace(Φ(t)^T·Φ(t)) / n
 *
 * This is the average squared gain from initial condition to state.
 * Large μ indicates high sensitivity to initial conditions.
 *
 * @param Phi state transition matrix
 * @param n dimension
 * @return trace-based sensitivity measure
 */
double initial_condition_sensitivity_measure(const double *Phi, int n) {
    if (Phi == NULL || n <= 0) return 0.0;

    double trace = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            trace += Phi[j * n + i] * Phi[j * n + i];
        }
    }
    return trace / (double)n;
}

/* ==========================================================================
 * L5: Discrete-Time Trajectory Sensitivity
 *
 * For discrete-time system x_{k+1} = f(x_k, p):
 *   s_{k+1} = ∂f/∂x|_{x_k} · s_k + ∂f/∂p|_{x_k}
 *   s₀ = ∂x₀/∂p
 *
 * This is the discrete variational equation.
 * ========================================================================== */

/**
 * Integrate discrete-time trajectory sensitivity.
 *
 * @param f_dt function: x_next = f(x, p)
 * @param n state dimension
 * @param n_params number of parameters
 * @param params parameter values
 * @param x0 initial state
 * @param dx0_dp initial state sensitivity (n×n_params) or NULL
 * @param n_steps number of time steps
 * @param x_traj output: state trajectory ((n_steps+1)×n, row-major)
 * @param s_traj output: sensitivity trajectory (n_params × (n_steps+1)×n)
 */
void discrete_trajectory_sensitivity(
    void (*f_dt)(const double *x, const double *p, int n, int np, double *x_next),
    int n, int n_params,
    const double *params, const double *x0, const double *dx0_dp,
    int n_steps, double *x_traj, double *s_traj) {

    if (f_dt == NULL || x_traj == NULL || s_traj == NULL || n <= 0) return;

    /* Initialize state */
    for (int i = 0; i < n; i++) x_traj[i] = x0[i];

    /* Initialize sensitivity */
    if (dx0_dp != NULL) {
        for (int p = 0; p < n_params; p++) {
            for (int i = 0; i < n; i++) {
                s_traj[p * (n_steps + 1) * n + i] = dx0_dp[p * n + i];
            }
        }
    } else {
        for (int p = 0; p < n_params; p++) {
            for (int i = 0; i < n; i++) {
                s_traj[p * (n_steps + 1) * n + i] = 0.0;
            }
        }
    }

    double *x_curr = (double *)malloc((size_t)n * sizeof(double));
    double *x_next = (double *)malloc((size_t)n * sizeof(double));
    double *Jx = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *Jp = (double *)malloc((size_t)(n * n_params) * sizeof(double));

    if (x_curr == NULL || x_next == NULL || Jx == NULL || Jp == NULL) {
        free(x_curr); free(x_next); free(Jx); free(Jp);
        return;
    }

    memcpy(x_curr, x0, (size_t)n * sizeof(double));

    for (int k = 0; k < n_steps; k++) {
        /* Compute x_{k+1} = f(x_k, p) */
        f_dt(x_curr, params, n, n_params, x_next);

        /* Estimate Jacobians via finite differences */
        double h = 1e-6;
        double *x_pert = (double *)malloc((size_t)n * sizeof(double));

        /* ∂f_i/∂x_j */
        for (int j = 0; j < n; j++) {
            memcpy(x_pert, x_curr, (size_t)n * sizeof(double));
            x_pert[j] += h;
            double *f_plus = (double *)malloc((size_t)n * sizeof(double));
            f_dt(x_pert, params, n, n_params, f_plus);

            x_pert[j] -= 2.0 * h;
            double *f_minus = (double *)malloc((size_t)n * sizeof(double));
            f_dt(x_pert, params, n, n_params, f_minus);

            for (int i = 0; i < n; i++) {
                Jx[i * n + j] = (f_plus[i] - f_minus[i]) / (2.0 * h);
            }

            free(f_plus); free(f_minus);
        }

        /* ∂f_i/∂p_j */
        double *p_pert = (double *)malloc((size_t)n_params * sizeof(double));
        for (int j = 0; j < n_params; j++) {
            memcpy(p_pert, params, (size_t)n_params * sizeof(double));
            p_pert[j] += h;
            double *f_p = (double *)malloc((size_t)n * sizeof(double));
            f_dt(x_curr, p_pert, n, n_params, f_p);

            p_pert[j] -= 2.0 * h;
            double *f_m = (double *)malloc((size_t)n * sizeof(double));
            f_dt(x_curr, p_pert, n, n_params, f_m);

            for (int i = 0; i < n; i++) {
                Jp[i * n_params + j] = (f_p[i] - f_m[i]) / (2.0 * h);
            }

            free(f_p); free(f_m);
        }

        /* Store state */
        for (int i = 0; i < n; i++) {
            x_traj[(k + 1) * n + i] = x_next[i];
        }

        /* Update sensitivities: s_{k+1} = Jx·s_k + Jp */
        for (int p = 0; p < n_params; p++) {
            double *s_k = &s_traj[p * (n_steps + 1) * n + k * n];
            double *s_kp1 = &s_traj[p * (n_steps + 1) * n + (k + 1) * n];

            for (int i = 0; i < n; i++) {
                s_kp1[i] = 0.0;
                for (int j = 0; j < n; j++) {
                    s_kp1[i] += Jx[i * n + j] * s_k[j];
                }
                s_kp1[i] += Jp[i * n_params + p];
            }
        }

        /* Advance state */
        memcpy(x_curr, x_next, (size_t)n * sizeof(double));

        free(x_pert);
        free(p_pert);
    }

    free(x_curr); free(x_next); free(Jx); free(Jp);
}

/* ==========================================================================
 * L8: Periodic Orbit Sensitivity via Monodromy Matrix
 *
 * For a periodic orbit with period T: x(t+T) = x(t),
 * the monodromy matrix M = Φ(T, 0) is the state transition matrix
 * evaluated over one period.
 *
 * Floquet multipliers μ_i are the eigenvalues of M.
 * Floquet exponents λ_i = ln(μ_i)/T.
 *
 * Sensitivity of Floquet multipliers to parameters:
 *   ∂μ_k/∂p = y_k^* · (∂M/∂p) · v_k / (y_k^*·v_k)
 * where v_k = right eigenvector, y_k = left eigenvector of M.
 *
 * This follows from Wilkinson's formula applied to the monodromy matrix.
 * ========================================================================== */

/**
 * Compute the monodromy matrix M = Φ(T, 0) for a periodic system.
 *
 * @param A_periodic function: A(t+T) = A(t)
 * @param period T
 * @param n dimension
 * @param n_steps integration steps per period
 * @param M output: monodromy matrix (n×n)
 */
void monodromy_matrix(
    void (*A_periodic)(double t, int n, double *A),
    double period, int n, int n_steps, double *M) {

    if (A_periodic == NULL || M == NULL || n <= 0) return;

    double dt = period / (double)n_steps;

    /* Initialize M = I */
    for (int i = 0; i < n * n; i++) M[i] = 0.0;
    for (int i = 0; i < n; i++) M[i * n + i] = 1.0;

    double *A = (double *)malloc((size_t)(n * n) * sizeof(double));
    double *M_new = (double *)malloc((size_t)(n * n) * sizeof(double));

    if (A == NULL || M_new == NULL) {
        free(A); free(M_new);
        return;
    }

    /* Euler integration: M_{k+1} = (I + dt·A(t_k)) · M_k */
    for (int step = 0; step < n_steps; step++) {
        double t = (double)step * dt;
        A_periodic(t, n, A);

        /* M_new = (I + dt·A)·M */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                M_new[i * n + j] = M[i * n + j];
                for (int k = 0; k < n; k++) {
                    M_new[i * n + j] += dt * A[i * n + k] * M[k * n + j];
                }
            }
        }

        memcpy(M, M_new, (size_t)(n * n) * sizeof(double));
    }

    free(A); free(M_new);
}

/**
 * Compute Floquet multipliers (eigenvalues of monodromy matrix).
 *
 * @param M monodromy matrix (n×n)
 * @param n dimension
 * @param multipliers output: Floquet multipliers (n Complex numbers)
 * @return 0 on success
 */
int floquet_multipliers(const double *M, int n, Complex *multipliers) {
    if (M == NULL || multipliers == NULL || n <= 0) return -1;

    double *M_copy = (double *)malloc((size_t)(n * n) * sizeof(double));
    if (M_copy == NULL) return -1;

    memcpy(M_copy, M, (size_t)(n * n) * sizeof(double));
    int ret = qr_eigenvalues(M_copy, n, multipliers);
    free(M_copy);
    return ret;
}

/**
 * Check orbital stability from Floquet multipliers.
 *
 * A periodic orbit is asymptotically stable if all Floquet multipliers
 * lie strictly inside the unit circle (|μ_i| < 1), except one trivial
 * multiplier at μ = 1 corresponding to the tangential direction.
 *
 * @param multipliers Floquet multipliers
 * @param n number of multipliers
 * @param tol tolerance for treating |μ-1| as the trivial multiplier
 * @return 1 if orbit is stable, 0 otherwise
 */
int is_orbit_stable(const Complex *multipliers, int n, double tol) {
    if (multipliers == NULL || n <= 0) return 0;

    int trivial_found = 0;

    for (int i = 0; i < n; i++) {
        double abs_val = complex_abs(multipliers[i]);

        /* Check for trivial multiplier μ ≈ 1 */
        Complex diff = complex_sub(multipliers[i], complex_make(1.0, 0.0));
        if (complex_abs(diff) < tol && !trivial_found) {
            trivial_found = 1;
            continue;
        }

        if (abs_val >= 1.0) return 0;
    }

    return 1;
}

/**
 * Sensitivity of Floquet multiplier to parameter using Wilkinson's formula.
 *
 * ∂μ_k/∂p = (y_k^* · (∂M/∂p) · v_k) / (y_k^*·v_k)
 *
 * @param M monodromy matrix (n×n)
 * @param dM_dp sensitivity of M to parameter (n×n)
 * @param n dimension
 * @param mu_k which Floquet multiplier (0-indexed)
 * @param multipliers all Floquet multipliers
 * @param right_vectors right eigenvectors of M (n×n, row-major)
 * @param left_vectors left eigenvectors of M (n×n, row-major)
 * @return ∂μ_k/∂p
 */
Complex floquet_sensitivity(const double *M, const double *dM_dp,
                            int n, int mu_k,
                            const Complex *multipliers,
                            const double *right_vectors,
                            const double *left_vectors) {
    (void)M;  /* M used only indirectly via multipliers */
    if (dM_dp == NULL || right_vectors == NULL || left_vectors == NULL) {
        return complex_make(NAN, NAN);
    }

    double *x = (double *)malloc((size_t)n * sizeof(double));
    double *y = (double *)malloc((size_t)n * sizeof(double));

    if (x == NULL || y == NULL) {
        free(x); free(y);
        return complex_make(NAN, NAN);
    }

    /* Extract eigenvectors (assume real for simplicity) */
    for (int i = 0; i < n; i++) {
        x[i] = right_vectors[mu_k * n + i];
        y[i] = left_vectors[mu_k * n + i];
    }

    Complex mu_k_val = multipliers[mu_k];
    Complex sens = eigenvalue_param_sensitivity(n, dM_dp, mu_k_val, x, y);

    free(x); free(y);
    return sens;
}
