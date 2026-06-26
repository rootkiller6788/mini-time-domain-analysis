#ifndef DISTURBANCE_REJECTION_H
#define DISTURBANCE_REJECTION_H

#include "steady_state_error.h"

/* ============================================================================
 * mini-steady-state-error: Disturbance Rejection Module
 *
 * L2: Disturbance steady-state error analysis
 * L4: Superposition principle for multiple inputs
 * L5: Disturbance observer design
 * L6: Canonical disturbance rejection problems
 *
 * Reference: Nise §7.6, Ogata §5.6, Franklin §4.4
 * ============================================================================ */

/* L2: Disturbance input types */
typedef enum {
    DIST_STEP = 0,        /* constant disturbance (bias, offset) */
    DIST_RAMP = 1,        /* drifting disturbance */
    DIST_SINUSOIDAL = 2,  /* periodic disturbance (vibration, ripple) */
    DIST_IMPULSE = 3      /* shock/impulse disturbance */
} DisturbanceType;

/* L2: Disturbance injection point */
typedef enum {
    DIST_INPUT = 0,       /* disturbance at plant input (load disturbance) */
    DIST_OUTPUT = 1,      /* disturbance at plant output (sensor noise) */
    DIST_PLANT = 2        /* disturbance inside plant dynamics */
} DisturbanceLocation;

/* L2: Disturbance steady-state error result */
typedef struct {
    double e_ss_dist;          /* steady-state error due to disturbance */
    double e_ss_total;         /* total e_ss (reference + disturbance) */
    double rejection_ratio;    /* |e_ss_with_dist| / |e_ss_without_dist| */
    double disturbance_dc_gain;/* DC gain from disturbance to error */
    bool   disturbance_rejected; /* true if e_ss_dist = 0 */
} DisturbanceError;

/* L2: Disturbance-to-output transfer function result */
typedef struct {
    TransferFunction tf_dist_to_output;  /* Y(s)/D(s) */
    TransferFunction tf_dist_to_error;   /* E(s)/D(s) = -H(s)*Y(s)/D(s) for unity FB */
    double dc_gain_from_dist;            /* lim_{s→0} E(s)/D(s) */
} DisturbanceTF;

/* ---- L2: Core disturbance error computation ---- */

/**
 * L2: Compute steady-state error due to a disturbance.
 *
 * For disturbance D(s) at plant input (unity feedback, H=1):
 *   E(s) = -G(s) / (1 + G(s)) * D(s)
 *   e_ss_dist = lim_{s→0} s * E(s)
 *
 * @param G         plant transfer function (forward path)
 * @param dist_type type of disturbance
 * @param location  where disturbance enters the loop
 * @param magnitude disturbance magnitude
 * @return          disturbance steady-state error
 *
 * Complexity: O(n) where n = denom_order of G(s)
 */
double drej_compute_disturbance_error(const TransferFunction *G,
                                       DisturbanceType dist_type,
                                       DisturbanceLocation location,
                                       double magnitude);

/**
 * L2: Compute disturbance-to-error transfer function.
 *
 * For disturbance at plant input (unity feedback):
 *   E(s)/D(s) = -G(s) / (1 + G(s))
 *
 * For disturbance at plant output (unity feedback):
 *   E(s)/D(s) = -1 / (1 + G(s))
 *
 * Complexity: O((n+m)^2) polynomial multiplication
 */
DisturbanceTF drej_disturbance_transfer_function(const TransferFunction *G,
                                                   DisturbanceLocation location);

/**
 * L2: Check if disturbance is rejected in steady state.
 * Disturbance rejection condition: the loop must have a model of the
 * disturbance generator (Internal Model Principle applied to disturbance).
 *
 * For step disturbance at plant input: need integrator BEFORE disturbance entry
 * For step disturbance at plant output: need integrator in forward path
 *
 * Complexity: O(n) where n = denom_order
 */
bool drej_is_disturbance_rejected(const TransferFunction *G,
                                   DisturbanceType dist_type,
                                   DisturbanceLocation location);

/**
 * L4: Superposition applied to simultaneous reference + disturbance.
 * Total error = error from reference + error from disturbance.
 *
 * e_ss_total = e_ss_ref + e_ss_dist
 *
 * Complexity: O(n)
 */
DisturbanceError drej_combined_error(const TransferFunction *G,
                                      TestInputType ref_type, double ref_magnitude,
                                      DisturbanceType dist_type,
                                      DisturbanceLocation dist_location,
                                      double dist_magnitude);

/* ---- L5: Disturbance observer and compensation ---- */

/**
 * L5: Design disturbance observer gain for input disturbance estimation.
 *
 * Disturbance observer (DOB) estimates the disturbance and compensates
 * through feedforward cancellation.
 *
 * Given plant G(s) = b0 / (a0 + a1*s + a2*s^2 + ...)  (assumed stable),
 * the Q-filter cutoff determines the trade-off between disturbance rejection
 * bandwidth and noise sensitivity.
 *
 * @param G               plant transfer function
 * @param Q_cutoff_freq   Q-filter cutoff frequency (rad/s)
 * @param observer_gain   output: observer gain
 *
 * Complexity: O(n)
 */
void drej_design_disturbance_observer(const TransferFunction *G,
                                       double Q_cutoff_freq,
                                       double *observer_gain);

/**
 * L5: Compute improvement in disturbance rejection after adding DOB.
 * Returns the ratio of e_ss_dist_without_DOB / e_ss_dist_with_DOB.
 *
 * Complexity: O(1)
 */
double drej_observer_improvement(double original_e_ss_dist,
                                  double observer_gain,
                                  double plant_dc_gain);

/**
 * L5: Feedforward disturbance compensation.
 * If disturbance is measurable, use feedforward to cancel its effect.
 *
 * Given disturbance-to-output TF G_d(s) and control-to-output TF G_u(s),
 * the ideal feedforward is: G_ff(s) = -G_d(s) / G_u(s)
 *
 * @param G_u     control-to-output transfer function
 * @param G_d     disturbance-to-output transfer function
 * @param G_ff    output: feedforward compensator (caller frees)
 *
 * Complexity: O((n+m)^2)
 */
void drej_feedforward_design(const TransferFunction *G_u,
                              const TransferFunction *G_d,
                              TransferFunction *G_ff);

/* ---- L6: Canonical disturbance rejection problems ---- */

/**
 * L6: DC motor load torque disturbance rejection.
 * For DC motor: G(s) = K / (s*(tau*s + 1))
 * Load torque T_L acts as input disturbance.
 * Computes steady-state speed error due to load torque.
 *
 * @param K_motor       motor gain (rad/s/V)
 * @param tau_motor     motor time constant (s)
 * @param T_load        load torque magnitude (N*m)
 * @return              steady-state speed error (rad/s)
 *
 * Complexity: O(1)
 */
double drej_dc_motor_load_torque(double K_motor, double tau_motor, double T_load);

/**
 * L6: Elevator load disturbance.
 * For elevator system with PI controller, compute steady-state position
 * error when a constant load (passenger weight) is applied.
 *
 * @param K_p           proportional gain
 * @param K_i           integral gain
 * @param load_mass     load mass (kg)
 * @param motor_const   motor force constant (N/V)
 * @return              steady-state position error (m)
 *
 * Complexity: O(1)
 */
double drej_elevator_load_error(double K_p, double K_i,
                                 double load_mass, double motor_const);

/**
 * L6: Temperature control disturbance rejection.
 * For thermal system: G(s) = K_th / (tau_th*s + 1)
 * Environmental temperature T_amb acts as output disturbance.
 * PI controller: G_c(s) = K_p + K_i/s
 *
 * @param K_th          thermal system gain
 * @param tau_th        thermal time constant (s)
 * @param K_p           PI proportional gain
 * @param K_i           PI integral gain
 * @param T_amb         ambient temperature disturbance (deg C)
 * @return              steady-state temperature error (deg C)
 *
 * Complexity: O(1)
 */
double drej_thermal_disturbance(double K_th, double tau_th,
                                 double K_p, double K_i, double T_amb);

/* ---- Utility ---- */

void drej_error_print(const DisturbanceError *err);
const char *drej_type_name(DisturbanceType t);
const char *drej_location_name(DisturbanceLocation loc);

#endif /* DISTURBANCE_REJECTION_H */
