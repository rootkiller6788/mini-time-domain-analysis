/* ============================================================================
 * mini-steady-state-error: Disturbance Rejection Implementation
 * L2: Disturbance steady-state error analysis
 * L4: Superposition principle for reference + disturbance
 * L5: Disturbance observer and feedforward design
 * L6: Canonical problems (DC motor, elevator, thermal)
 * Reference: Nise 7.6, Ogata 5.6, Franklin 4.4
 * ============================================================================ */

#include "disturbance_rejection.h"
#include "error_constants.h"
#include "system_type.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

double drej_compute_disturbance_error(const TransferFunction *G,
                                       DisturbanceType dist_type,
                                       DisturbanceLocation location,
                                       double magnitude)
{
    if (!G) return NAN;
    double G0 = sse_dc_gain(G);
    bool has_int = (sse_determine_system_type(G) >= SYSTEM_TYPE_1);
    double e_ss = 0.0;

    if (location == DIST_INPUT) {
        switch (dist_type) {
        case DIST_STEP:
            if (has_int) e_ss = 0.0;
            else if (fabs(1.0+G0)<1e-15) e_ss = -INFINITY;
            else e_ss = -G0*magnitude/(1.0+G0);
            break;
        case DIST_RAMP:
            if (sse_determine_system_type(G)>=SYSTEM_TYPE_2) e_ss=0.0;
            else if (has_int) {
                double Kv=ec_compute_Kv(G);
                if (fabs(Kv)>1e-15) e_ss=-magnitude/Kv;
                else e_ss=-INFINITY;
            } else e_ss=-INFINITY;
            break;
        case DIST_SINUSOIDAL: e_ss=0.0; break;
        case DIST_IMPULSE: e_ss=0.0; break;
        }
    } else if (location == DIST_OUTPUT) {
        switch (dist_type) {
        case DIST_STEP:
            if (has_int) e_ss=0.0;
            else e_ss=-magnitude/(1.0+G0);
            break;
        case DIST_RAMP:
            if (sse_determine_system_type(G)>=SYSTEM_TYPE_2) e_ss=0.0;
            else if (has_int) {
                double Kv=ec_compute_Kv(G);
                if (fabs(Kv)>1e-15) e_ss=-magnitude/(Kv*G0);
                else e_ss=-INFINITY;
            } else e_ss=-INFINITY;
            break;
        case DIST_SINUSOIDAL: e_ss=0.0; break;
        case DIST_IMPULSE: e_ss=0.0; break;
        }
    } else {
        e_ss = magnitude;
    }
    return e_ss;
}

DisturbanceTF drej_disturbance_transfer_function(const TransferFunction *G,
                                                   DisturbanceLocation location)
{
    DisturbanceTF result;
    memset(&result, 0, sizeof(result));
    if (!G) return result;

    TransferFunction one_plus_G = tf_alloc(G->numer_order,
        (G->numer_order > G->denom_order) ? G->numer_order : G->denom_order);
    if (one_plus_G.numer_order < 0) return result;

    one_plus_G.gain = 1.0;
    int i;
    for (i = 0; i <= G->denom_order; i++)
        one_plus_G.denom_coeffs[i] = G->denom_coeffs[i];
    for (i = 0; i <= G->numer_order; i++)
        one_plus_G.numer_coeffs[i] = G->denom_coeffs[i]
            + G->gain * G->numer_coeffs[i];

    if (location == DIST_INPUT) {
        result.tf_dist_to_error = tf_alloc(G->numer_order, one_plus_G.numer_order);
        if (result.tf_dist_to_error.numer_order >= 0) {
            result.tf_dist_to_error.gain = -G->gain;
            for (i = 0; i <= G->numer_order; i++)
                result.tf_dist_to_error.numer_coeffs[i] = G->numer_coeffs[i];
            for (i = 0; i <= one_plus_G.numer_order; i++)
                result.tf_dist_to_error.denom_coeffs[i] = one_plus_G.numer_coeffs[i];
        }
    } else {
        result.tf_dist_to_error = tf_alloc(0, one_plus_G.numer_order);
        if (result.tf_dist_to_error.numer_order >= 0) {
            result.tf_dist_to_error.gain = -1.0;
            result.tf_dist_to_error.numer_coeffs[0] = 1.0;
            for (i = 0; i <= one_plus_G.numer_order; i++)
                result.tf_dist_to_error.denom_coeffs[i] = one_plus_G.numer_coeffs[i];
        }
    }

    result.dc_gain_from_dist = sse_dc_gain(&result.tf_dist_to_error);
    tf_free(&one_plus_G);
    return result;
}

bool drej_is_disturbance_rejected(const TransferFunction *G,
                                   DisturbanceType dist_type,
                                   DisturbanceLocation location)
{
    if (!G) return false;
    int sys_type = (int)sse_determine_system_type(G);
    (void)location;
    switch (dist_type) {
    case DIST_STEP:   return (sys_type >= 1);
    case DIST_RAMP:   return (sys_type >= 2);
    case DIST_SINUSOIDAL: return false;
    case DIST_IMPULSE:    return true;
    default: return false;
    }
}

DisturbanceError drej_combined_error(const TransferFunction *G,
                                      TestInputType ref_type, double ref_mag,
                                      DisturbanceType dist_type,
                                      DisturbanceLocation dist_loc,
                                      double dist_mag)
{
    DisturbanceError result;
    memset(&result, 0, sizeof(result));
    if (!G) return result;

    double Kp = ec_compute_Kp(G), Kv = ec_compute_Kv(G), Ka = ec_compute_Ka(G);
    double e_ss_ref = sse_compute_specific(sse_determine_system_type(G),
        ref_type, Kp, Kv, Ka, ref_mag);
    double e_ss_dist = drej_compute_disturbance_error(G, dist_type, dist_loc, dist_mag);

    result.e_ss_dist = e_ss_dist;
    result.e_ss_total = e_ss_ref + e_ss_dist;
    result.disturbance_rejected = (fabs(e_ss_dist) < 1e-10);
    if (fabs(e_ss_ref) > 1e-10)
        result.rejection_ratio = fabs(e_ss_dist / e_ss_ref);
    else
        result.rejection_ratio = (fabs(e_ss_dist) < 1e-10) ? 1.0 : INFINITY;
    result.disturbance_dc_gain = sse_dc_gain(G) / (1.0 + sse_dc_gain(G));
    return result;
}

void drej_design_disturbance_observer(const TransferFunction *G,
                                       double Q_cutoff_freq,
                                       double *observer_gain)
{
    if (!G || !observer_gain) return;
    double G0 = sse_dc_gain(G);
    if (fabs(G0) < 1e-15) { *observer_gain = Q_cutoff_freq; return; }
    *observer_gain = Q_cutoff_freq / G0;
}

double drej_observer_improvement(double original_e_ss_dist,
                                  double observer_gain, double plant_dc_gain)
{
    if (fabs(original_e_ss_dist) < 1e-15) return 1.0;
    return 1.0 + observer_gain * fabs(plant_dc_gain);
}

void drej_feedforward_design(const TransferFunction *G_u,
                              const TransferFunction *G_d,
                              TransferFunction *G_ff)
{
    if (!G_u || !G_d || !G_ff) return;
    double G_u0 = sse_dc_gain(G_u);
    double G_d0 = sse_dc_gain(G_d);
    if (fabs(G_u0) < 1e-15) {
        *G_ff = tf_alloc(0, 0);
        if (G_ff->numer_order >= 0) G_ff->numer_coeffs[0] = 0.0;
        return;
    }
    *G_ff = tf_alloc(0, 0);
    if (G_ff->numer_order >= 0) {
        G_ff->gain = 1.0; G_ff->numer_coeffs[0] = -G_d0 / G_u0;
        G_ff->denom_coeffs[0] = 1.0;
    }
}

/* ============================================================================
 * L6: Canonical Disturbance Rejection Problems
 * ============================================================================ */

double drej_dc_motor_load_torque(double K_motor, double tau_motor, double T_load)
{
    /* DC motor: G(s) = K_motor/(s*(tau_motor*s+1))
     * Load torque T_L acts as input disturbance.
     * With unity FB: E(s)/T_L(s) = -G(s)/(1+G(s))
     * For step T_L: e_ss = lim s*(-G(s)*(T_L/s)/(1+G(s))) = -G(0)*T_L/(1+G(0))
     * G(0) = inf (Type 1) so e_ss = 0.
     * With PI controller, the integrator rejects constant load torque.
     * Without integrator (Type 0 with P control):
     *   Steady-state speed error = T_L * R_a / (K_t * K_b + R_a * B)
     *   where R_a=armature resistance, K_t=torque const, K_b=back EMF const.
     * This function returns the speed error under P control with gain K_motor. */
    if (fabs(K_motor) < 1e-15) return INFINITY;
    /* For Type 1 system (motor has integrator), constant disturbance is rejected */
    (void)tau_motor; (void)T_load;
    return 0.0; /* motor has inherent integrator */
}

double drej_elevator_load_error(double K_p, double K_i,
                                 double load_mass, double motor_const)
{
    /* Elevator: PI-controlled position system.
     * Load = passenger weight acts as output disturbance.
     * G_plant(s) = motor_const / (M*s^2 + B*s)  [simplified]
     * PI: G_c(s) = K_p + K_i/s
     * For step load disturbance: e_ss = -load/G_ol(0) where G_ol(0) depends on
     * the type of the loop.
     * With PI (Type 1+), constant output disturbance is rejected: e_ss = 0. */
    if (fabs(K_i) < 1e-15) {
        /* P-only: e_ss = load_mass * g / (motor_const * K_p) */
        double g = 9.81;
        return load_mass * g / (motor_const * K_p);
    }
    /* With integral action: zero steady-state error */
    (void)load_mass; (void)motor_const;
    return 0.0;
}

double drej_thermal_disturbance(double K_th, double tau_th,
                                 double K_p, double K_i, double T_amb)
{
    /* Thermal system: G(s) = K_th/(tau_th*s+1)
     * PI controller: G_c(s) = K_p + K_i/s
     * Ambient temperature T_amb acts as output disturbance.
     * e_ss_dist = lim s*[-T_amb/s * 1/(1+G_c(s)G(s))]
     *           = -T_amb / (1 + G_c(0)*G(0))
     * With PI: G_c(0) = inf (integrator), so e_ss_dist = 0.
     * With P-only: G_c(0) = K_p, e_ss_dist = -T_amb/(1+K_p*K_th) */
    if (fabs(K_i) < 1e-15) {
        double G0 = K_p * K_th;
        return -T_amb / (1.0 + G0);
    }
    (void)tau_th; (void)K_p;
    return 0.0; /* PI rejects constant output disturbance */
}

void drej_error_print(const DisturbanceError *err)
{
    if (!err) return;
    printf("Disturbance Error Analysis:\n");
    printf("  e_ss (disturbance): %.6g\n", err->e_ss_dist);
    printf("  e_ss (total): %.6g\n", err->e_ss_total);
    printf("  Rejection ratio: %.6g\n", err->rejection_ratio);
    printf("  Disturbance rejected: %s\n",
           err->disturbance_rejected ? "YES" : "NO");
}

const char *drej_type_name(DisturbanceType t) {
    switch(t){case DIST_STEP:return "Step";case DIST_RAMP:return "Ramp";
    case DIST_SINUSOIDAL:return "Sinusoidal";case DIST_IMPULSE:return "Impulse";
    default:return "Unknown";}
}
const char *drej_location_name(DisturbanceLocation loc) {
    switch(loc){case DIST_INPUT:return "Plant Input";case DIST_OUTPUT:return "Plant Output";
    case DIST_PLANT:return "Internal Plant";default:return "Unknown";}
}
