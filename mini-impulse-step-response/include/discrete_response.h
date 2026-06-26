/**
 * @file discrete_response.h
 * @brief Discrete-time impulse and step response computation.
 *
 * For digital control systems and sampled-data systems, the discrete-time
 * impulse response h[k] and step response s[k] characterize the system
 * behavior in the Z-domain.
 *
 * The discrete unit impulse (Kronecker delta) is:
 *   delta[k] = 1 for k=0, 0 otherwise
 *
 * The discrete unit step is:
 *   u[k] = 1 for k >= 0, 0 otherwise
 *
 * L1 --- Definitions: Discrete impulse delta[k], discrete step u[k],
 *       discrete impulse response h[k], discrete step response s[k]
 * L2 --- Core Concepts: Z-transform relationship, discrete convolution,
 *       relationship between continuous and discrete time
 * L4 --- Fundamental Laws: Z{h[k]} = H(z), discrete convolution sum,
 *       s[k] = sum_{i=0}^{k} h[i]
 * L5 --- Engineering Methods: Difference equation simulation,
 *       zero-order hold (ZOH) discretization
 *
 * Course alignment:
 *   MIT 6.302 -- digital control, discrete-time response
 *   Stanford ENGR105 -- sampled-data systems
 *   Berkeley ME132 -- discrete-time dynamic systems
 *   ETH 151-0591 -- zeitdiskrete Systeme
 *   Tsinghua -- discrete system analysis
 */

#ifndef DISCRETE_RESPONSE_H
#define DISCRETE_RESPONSE_H

#include "time_response_common.h"

/* ==========================================================================
 * L1: Discrete-Time System Models
 * ========================================================================== */

/**
 * Discrete-time transfer function (pulse transfer function):
 *   H(z) = (b0 + b1*z^{-1} + ... + bm*z^{-m}) / (1 + a1*z^{-1} + ... + an*z^{-n})
 *
 * This is the Z-transform of the discrete impulse response h[k].
 */
typedef struct {
    int     num_deg;  /**< numerator degree m */
    int     den_deg;  /**< denominator degree n */
    double *num;      /**< numerator coefficients [b0, b1, ..., bm] */
    double *den;      /**< denominator coefficients [1, a1, ..., an] (den[0]=1) */
} DiscreteTransferFunction;

/**
 * Discrete state-space model:
 *   x[k+1] = Ad * x[k] + Bd * u[k]
 *   y[k]   = Cd * x[k] + Dd * u[k]
 */
typedef struct {
    int     n_states;
    int     n_inputs;
    int     n_outputs;
    double *Ad;  /**< discrete system matrix [n*n] */
    double *Bd;  /**< discrete input matrix [n*m] */
    double *Cd;  /**< discrete output matrix [p*n] */
    double *Dd;  /**< discrete feedthrough [p*m] */
} DiscreteStateSpace;

/* ==========================================================================
 * L1: Discrete Response Computation
 * ========================================================================== */

/**
 * Compute discrete impulse response h[k] for a discrete transfer function.
 *
 * Method: long division of numerator by denominator of H(z).
 * The coefficients of the series H(z) = h[0] + h[1]*z^{-1} + h[2]*z^{-2} + ...
 * are the impulse response values.
 *
 * @param dtf   Discrete transfer function.
 * @param N     Number of samples to compute.
 * @param h_out Output: impulse response array (size N, caller-allocated).
 *
 * Complexity: O(N * den_deg)
 * Textbook: Ogata Section 4-4
 */
void discrete_impulse_response(const DiscreteTransferFunction *dtf,
                                int N, double *h_out);

/**
 * Compute discrete step response s[k] for a discrete transfer function.
 *
 * Method 1 (direct): simulate difference equation with u[k] = 1.
 * Method 2 (via h[k]): s[k] = sum_{i=0}^{k} h[i] (cumulative sum).
 *
 * This function uses Method 1 (difference equation) for numerical stability.
 *
 * @param dtf   Discrete transfer function.
 * @param N     Number of samples.
 * @param s_out Output: step response array (size N).
 *
 * Complexity: O(N * den_deg)
 */
void discrete_step_response(const DiscreteTransferFunction *dtf,
                             int N, double *s_out);

/**
 * Simulate a discrete state-space model with an arbitrary input sequence.
 *
 * x[0] = 0, then:
 *   x[k+1] = Ad*x[k] + Bd*u[k]
 *   y[k]   = Cd*x[k] + Dd*u[k]
 *
 * @param dss   Discrete state-space model.
 * @param u     Input sequence (N values).
 * @param N     Number of time steps.
 * @return ResponseTrajectory with y[k] values at times k*Ts.
 */
ResponseTrajectory *discrete_state_space_simulate(const DiscreteStateSpace *dss,
                                                    const double *u, int N);

/* ==========================================================================
 * L5: Continuous-to-Discrete Conversion
 * ========================================================================== */

/**
 * Discretize a continuous transfer function using Zero-Order Hold (ZOH).
 *
 * G(z) = (1 - z^{-1}) * Z{L^{-1}{G(s)/s} evaluated at t=k*Ts}
 *
 * For a first-order system G(s) = K/(tau*s + 1):
 *   G(z) = K * (1 - exp(-Ts/tau)) / (z - exp(-Ts/tau))
 *
 * @param tf   Continuous transfer function.
 * @param Ts   Sampling period [s].
 * @return Discrete transfer function (caller must discrete_tf_free).
 *
 * Complexity: O(n) for first/second order, O(n^3) for higher order
 * Textbook: Franklin, Powell, Workman "Digital Control" (1998)
 */
DiscreteTransferFunction *c2d_zoh(const TransferFunction *tf, double Ts);

/**
 * Discretize a continuous state-space model using ZOH.
 *
 * Ad = exp(A * Ts)
 * Bd = integral_0^{Ts} exp(A*(Ts-tau)) * B dtau
 *    = A^{-1} * (exp(A*Ts) - I) * B
 *
 * Uses matrix exponential (scaling-and-squaring + Pade) for Ad,
 * and numerical integration for Bd.
 *
 * @param ss  Continuous state-space model.
 * @param Ts  Sampling period [s].
 * @return Discrete state-space model.
 */
DiscreteStateSpace *c2d_ss_zoh(const StateSpaceModel *ss, double Ts);

/**
 * Discretize using Tustin's method (bilinear transform).
 *
 * s = (2/Ts) * (z-1)/(z+1)
 *
 * @param tf  Continuous transfer function.
 * @param Ts  Sampling period.
 * @return Discrete transfer function.
 */
DiscreteTransferFunction *c2d_tustin(const TransferFunction *tf, double Ts);

/* ==========================================================================
 * L2: Discrete Response Properties
 * ========================================================================== */

/**
 * Compute the final value of a discrete step response using
 * the discrete Final Value Theorem:
 *   lim_{k->inf} s[k] = lim_{z->1} (z-1)/z * H(z) * z/(z-1) = H(1)
 *
 * @param dtf  Discrete transfer function.
 * @return Steady-state step response value.
 */
double discrete_final_value(const DiscreteTransferFunction *dtf);

/**
 * Check stability of discrete system: all poles must satisfy |z| < 1.
 * Uses Jury's stability test for polynomial root location.
 *
 * @param dtf  Discrete transfer function.
 * @return 1 if stable, 0 otherwise.
 */
int discrete_is_stable(const DiscreteTransferFunction *dtf);

/**
 * Compute the discrete impulse response energy:
 * E_h = sum_{k=0}^{inf} (h[k])^2
 *
 * For finite N, approximates with first N terms.
 *
 * @param dtf  Discrete transfer function.
 * @param N    Number of terms to sum.
 * @return Approximate impulse energy.
 */
double discrete_impulse_energy(const DiscreteTransferFunction *dtf, int N);

/**
 * Compute settling time in discrete steps: smallest k such that
 * |s[k] - s[inf]| <= band * |s[inf]| for all m >= k.
 *
 * @param dtf   Discrete transfer function.
 * @param N     Maximum steps to check.
 * @param band  Settling band (e.g., 0.02).
 * @return Settling step index, or N if not settled.
 */
int discrete_settling_steps(const DiscreteTransferFunction *dtf,
                             int N, double band);

/* ==========================================================================
 * L2: Sampling Rate Selection
 * ========================================================================== */

/**
 * Recommend sampling frequency based on system bandwidth.
 * Rule: fs >= 20 * f_BW to fs >= 40 * f_BW for control applications.
 *
 * For second-order systems, the bandwidth is approximately:
 *   f_BW = (omega_n/(2*pi)) * sqrt(1 - 2*zeta^2 + sqrt(4*zeta^4 - 4*zeta^2 + 2))
 *
 * @param sys  Continuous second-order model.
 * @return Recommended sampling period Ts [s].
 */
double recommend_sampling_period(const SecondOrderModel *sys);

/**
 * Compute the aliasing error when sampling a continuous signal.
 * Compares the frequency response of continuous vs. discrete systems
 * at the Nyquist frequency.
 *
 * @param tf  Continuous transfer function.
 * @param Ts  Sampling period.
 * @return Maximum relative error in frequency response up to Nyquist.
 */
double aliasing_error_estimate(const TransferFunction *tf, double Ts);

/* ==========================================================================
 * L3: Memory Management
 * ========================================================================== */

void discrete_tf_free(DiscreteTransferFunction *dtf);
void discrete_ss_free(DiscreteStateSpace *dss);
DiscreteTransferFunction *discrete_tf_clone(const DiscreteTransferFunction *dtf);

#endif /* DISCRETE_RESPONSE_H */
