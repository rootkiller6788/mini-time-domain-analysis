/**
 * @file convolution_response.h
 * @brief Convolution integral computation for time-domain response analysis.
 *
 * The convolution integral is the fundamental operation for computing the
 * output of an LTI system given its impulse response h(t) and input u(t):
 *
 *   y(t) = (h * u)(t) = integral_0^t h(tau) * u(t - tau) dtau
 *                      = integral_0^t h(t - tau) * u(tau) dtau
 *
 * This is the time-domain equivalent of multiplication in the Laplace domain:
 *   Y(s) = H(s) * U(s)  <->  y(t) = h(t) * u(t)
 *
 * L1 --- Definitions: Convolution integral, commutative property,
 *       associative property, distributive property
 * L2 --- Core Concepts: LTI system input-output relationship,
 *       Green's function interpretation of h(t)
 * L4 --- Fundamental Laws: Convolution theorem L{f*g} = F(s)*G(s),
 *       BIBO stability criterion via h(t)
 * L5 --- Engineering Methods: Numerical quadrature, FFT-based convolution,
 *       recursive convolution for real-time implementation
 *
 * Course alignment:
 *   MIT 6.003 Signal Processing -- convolution
 *   Stanford ENGR105 -- LTI systems and convolution
 *   Berkeley ME132 -- convolution integral in dynamics
 *   Tsinghua -- convolution method
 */

#ifndef CONVOLUTION_RESPONSE_H
#define CONVOLUTION_RESPONSE_H

#include "time_response_common.h"

/* ==========================================================================
 * L5: Numerical Convolution Methods
 * ========================================================================== */

/**
 * Compute system output via direct numerical convolution of h(t) and u(t).
 *
 * y[k] = sum_{i=0}^{k} h[i] * u[k-i] * dt     (discrete convolution)
 *
 * This is the most general method: works for any h(t) and u(t) given as
 * sampled arrays. O(N^2) complexity where N is number of samples.
 *
 * @param h        Impulse response samples [N points].
 * @param u        Input signal samples [N points].
 * @param N        Number of samples.
 * @param dt       Sampling interval.
 * @return Allocated ResponseTrajectory (caller must response_trajectory_free).
 *
 * Complexity: O(N^2)
 * Textbook: Oppenheim & Willsky, "Signals and Systems" (1997)
 */
ResponseTrajectory *direct_convolution(const double *h, const double *u,
                                        int N, double dt);

/**
 * Compute convolution using the trapezoidal quadrature rule.
 *
 * y(t_k) = integral_0^{t_k} h(tau) * u(t_k - tau) dtau
 *
 * Each integral is evaluated via composite trapezoidal rule for O(h^2) accuracy.
 * Suitable for smooth h(t) and u(t).
 *
 * @param h        Impulse response samples.
 * @param u        Input signal samples.
 * @param N        Number of samples.
 * @param dt       Sampling interval.
 * @return ResponseTrajectory.
 *
 * Complexity: O(N^2)
 */
ResponseTrajectory *trapezoidal_convolution(const double *h, const double *u,
                                              int N, double dt);

/**
 * Compute convolution using Simpson's 1/3 rule.
 * Requires N to be odd (even number of intervals).
 * O(h^4) accuracy for smooth integrands.
 *
 * @param h        Impulse response samples (N odd).
 * @param u        Input signal samples (N odd).
 * @param N        Number of samples (must be odd).
 * @param dt       Sampling interval.
 * @return ResponseTrajectory, or NULL if N is even.
 *
 * Complexity: O(N^2)
 */
ResponseTrajectory *simpson_convolution(const double *h, const double *u,
                                         int N, double dt);

/**
 * Select the best quadrature method for convolution based on
 * signal characteristics (sample count, smoothness estimate).
 *
 * @param N  Number of samples.
 * @return Recommended QuadMethod.
 */
QuadMethod convolution_select_method(int N);

/* ==========================================================================
 * L5: Convolution with Analytic Impulse Response
 * ========================================================================== */

/**
 * Compute output by convolving an analytic first-order impulse response
 * with sampled input data. Uses the closed-form impulse response to avoid
 * sampling error in h(t).
 *
 * @param sys  First-order model.
 * @param u    Input signal samples.
 * @param N    Number of input samples.
 * @param dt   Sampling interval.
 * @return ResponseTrajectory.
 *
 * Complexity: O(N^2) with fast exp() evaluation
 */
ResponseTrajectory *first_order_convolution(const FirstOrderModel *sys,
                                              const double *u, int N, double dt);

/**
 * Compute output by convolving an analytic second-order impulse response
 * with sampled input data.
 *
 * Handles all damping cases (underdamped, critically damped, overdamped)
 * using the appropriate closed-form expression.
 *
 * @param sys  Second-order model.
 * @param u    Input signal samples.
 * @param N    Number of input samples.
 * @param dt   Sampling interval.
 * @return ResponseTrajectory.
 *
 * Complexity: O(N^2)
 */
ResponseTrajectory *second_order_convolution(const SecondOrderModel *sys,
                                               const double *u, int N, double dt);

/* ==========================================================================
 * L5: Fast Convolution Methods
 * ========================================================================== */

/**
 * Compute convolution using the overlap-add method with FFT.
 * For long signals, this is O(N log N) instead of O(N^2).
 *
 * Implements block-based FFT convolution:
 *   1. Partition input into blocks
 *   2. Zero-pad each block and h(t) to power-of-2 length
 *   3. FFT -> multiply -> IFFT
 *   4. Overlap-add to assemble output
 *
 * @param h         Impulse response samples.
 * @param h_len     Length of h.
 * @param u         Input signal samples.
 * @param u_len     Length of u.
 * @param dt        Sampling interval.
 * @return ResponseTrajectory.
 *
 * Complexity: O((h_len+u_len) * log(h_len+u_len))
 * Textbook: Smith, "Digital Signal Processing" (2003), Chapter 18
 */
ResponseTrajectory *fft_convolution(const double *h, int h_len,
                                     const double *u, int u_len, double dt);

/* ==========================================================================
 * L2: Convolution Properties and Verification
 * ========================================================================== */

/**
 * Verify the commutative property of convolution:
 * (h * u)(t) should equal (u * h)(t) to within numerical precision.
 *
 * @param h   Impulse response samples.
 * @param u   Input samples.
 * @param N   Number of samples.
 * @param dt  Sampling interval.
 * @param tol Tolerance for comparison.
 * @return 1 if commutative property holds, 0 otherwise.
 *
 * Complexity: O(N^2)
 */
int convolution_commutative_check(const double *h, const double *u,
                                   int N, double dt, double tol);

/**
 * Verify that convolution of h(t) with unit step u(t) equals
 * the step response s(t).
 *
 * Let u_step[k] = 1 for all k.
 * Then y_conv = h * u_step should equal step response s(t).
 *
 * @param h             Impulse response samples.
 * @param s_theoretical Theoretical step response samples for comparison.
 * @param N             Number of samples.
 * @param dt            Sampling interval.
 * @param tol           Tolerance.
 * @return Maximum absolute error.
 */
double convolution_step_verification(const double *h, const double *s_theoretical,
                                      int N, double dt);

/**
 * Verify that the area under h(t) equals DC gain via convolution with
 * a very long rectangular pulse approximating a step.
 *
 * @param h     Impulse response samples.
 * @param N     Number of samples.
 * @param dt    Sampling interval.
 * @param K_dc  Known DC gain.
 * @return Relative error |(integral h) - K_dc| / |K_dc|.
 */
double convolution_dc_gain_check(const double *h, int N, double dt, double K_dc);

/* ==========================================================================
 * L7: Convolution Applications
 * ========================================================================== */

/**
 * Compute the response of a structure to earthquake acceleration input
 * using convolution. The structure is modeled as a second-order SDOF
 * (single degree of freedom) system.
 *
 * @param sys         Second-order building model.
 * @param ground_acc  Ground acceleration time history [m/s^2].
 * @param N           Number of acceleration samples.
 * @param dt          Sampling interval.
 * @return Relative displacement response trajectory.
 */
ResponseTrajectory *seismic_response_convolution(const SecondOrderModel *sys,
                                                   const double *ground_acc,
                                                   int N, double dt);

/**
 * Compute the audio reverberation response using convolution.
 * The impulse response h(t) represents the room acoustics,
 * and the input u(t) is the dry audio signal.
 *
 * @param room_ir  Room impulse response (measured or synthetic).
 * @param ir_len   Length of room IR.
 * @param audio    Dry audio samples.
 * @param audio_len Length of audio.
 * @param dt       Sampling interval (1/sample_rate).
 * @return Reverberated audio response.
 */
ResponseTrajectory *reverberation_convolution(const double *room_ir, int ir_len,
                                                const double *audio, int audio_len,
                                                double dt);

#endif /* CONVOLUTION_RESPONSE_H */
