/**
 * @example ex3_vehicle_suspension.c
 * @brief Quarter-car suspension analysis: ride comfort and handling
 *
 * L6 --- Canonical Problem: Design and analyze a vehicle suspension
 * system using second-order models for both body bounce and wheel hop.
 * Evaluate ride comfort against ISO 2631 standards.
 *
 * L7 --- Application: Automotive suspension design, comparing
 * comfort-oriented (soft) vs. handling-oriented (stiff) tuning.
 *
 * Key trade-off: soft suspension → good comfort, poor handling
 *                stiff suspension → poor comfort, good handling
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/second_order.h"
#include "../include/transient_specs.h"
#include "../include/response_computation.h"
#include "../include/canonical_models.h"

int main(void)
{
    printf("=== Quarter-Car Suspension Analysis ===\n\n");

    /* Typical compact car parameters (Toyota Corolla-class) */
    QuarterCar comfort_tune = {
        .sprung_mass           = 350.0,   /* kg (corner) */
        .unsprung_mass         = 40.0,    /* kg (wheel assembly) */
        .suspension_stiffness  = 25000.0, /* N/m */
        .tire_stiffness        = 200000.0,/* N/m */
        .suspension_damping    = 1500.0   /* N·s/m */
    };

    /* Stiffer sports suspension */
    QuarterCar sport_tune = {
        .sprung_mass           = 350.0,
        .unsprung_mass         = 35.0,
        .suspension_stiffness  = 45000.0,
        .tire_stiffness        = 220000.0,
        .suspension_damping    = 2500.0
    };

    printf("Vehicle Configurations:\n");
    printf("%-20s  %12s  %12s\n", "Parameter", "Comfort", "Sport");
    printf("%-20s  %12s  %12s\n", "---------", "-------", "-----");
    printf("%-20s  %10.0f kg  %10.0f kg\n", "Sprung mass",
           comfort_tune.sprung_mass, sport_tune.sprung_mass);
    printf("%-20s  %10.0f kg  %10.0f kg\n", "Unsprung mass",
           comfort_tune.unsprung_mass, sport_tune.unsprung_mass);
    printf("%-20s  %10.0f N/m  %10.0f N/m\n", "Spring stiffness",
           comfort_tune.suspension_stiffness, sport_tune.suspension_stiffness);
    printf("%-20s  %10.0f N/m  %10.0f N/m\n", "Tire stiffness",
           comfort_tune.tire_stiffness, sport_tune.tire_stiffness);
    printf("%-20s  %10.0f Ns/m %10.0f Ns/m\n", "Damper",
           comfort_tune.suspension_damping, sport_tune.suspension_damping);

    /* Body bounce analysis */
    printf("\n=== Body Bounce Mode (Sprung Mass) ===\n");
    SecondOrderSystem body_comfort = quarter_car_body_model(&comfort_tune);
    SecondOrderSystem body_sport = quarter_car_body_model(&sport_tune);

    printf("%-25s  %12s  %12s\n", "Specification", "Comfort", "Sport");
    printf("%-25s  %12s  %12s\n", "------------", "-------", "-----");

    TransientSpecs spec_c = transient_compute_all(&body_comfort);
    TransientSpecs spec_s = transient_compute_all(&body_sport);

    printf("%-25s  %9.2f Hz  %9.2f Hz\n", "Natural freq (ωₙ/2π)",
           body_comfort.omega_n / (2.0 * M_PI),
           body_sport.omega_n / (2.0 * M_PI));
    printf("%-25s  %12.3f  %12.3f\n", "Damping ratio ζ",
           body_comfort.zeta, body_sport.zeta);
    printf("%-25s  %9.2f %%   %9.2f %%\n", "Overshoot PO",
           spec_c.percent_overshoot, spec_s.percent_overshoot);
    printf("%-25s  %9.3f s   %9.3f s\n", "Settling time (2%)",
           spec_c.settling_time_2pct, spec_s.settling_time_2pct);

    /* Wheel hop analysis */
    printf("\n=== Wheel Hop Mode (Unsprung Mass) ===\n");
    SecondOrderSystem wheel_comfort = quarter_car_wheel_model(&comfort_tune);
    SecondOrderSystem wheel_sport = quarter_car_wheel_model(&sport_tune);

    printf("%-25s  %12s  %12s\n", "Specification", "Comfort", "Sport");
    printf("%-25s  %12s  %12s\n", "------------", "-------", "-----");
    printf("%-25s  %9.2f Hz  %9.2f Hz\n", "Natural freq",
           wheel_comfort.omega_n / (2.0 * M_PI),
           wheel_sport.omega_n / (2.0 * M_PI));
    printf("%-25s  %12.3f  %12.3f\n", "Damping ratio ζ",
           wheel_comfort.zeta, wheel_sport.zeta);

    /* Ride comfort assessment (ISO 2631) */
    printf("\n=== Ride Comfort Assessment (ISO 2631 road profile) ===\n");
    printf("%-15s  %12s  %12s  %12s\n",
           "Freq [Hz]", "Comfort RMS", "Sport RMS", "ISO Limit");
    printf("%-15s  %12s  %12s  %12s\n",
           "---------", "-----------", "---------", "---------");

    double road_amplitude = 0.01;  /* 1 cm road roughness */
    double iso_limits[] = {0.315, 0.315, 0.5, 0.5, 0.8, 0.8};
    double freqs[] = {1.0, 2.0, 4.0, 8.0, 12.0, 16.0};

    for (int i = 0; i < 6; i++) {
        double rms_c = quarter_car_ride_comfort(
            &comfort_tune, freqs[i], road_amplitude);
        double rms_s = quarter_car_ride_comfort(
            &sport_tune, freqs[i], road_amplitude);
        printf("%-15.1f  %9.3f m/s² %9.3f m/s²  %9.3f m/s²",
               freqs[i], rms_c, rms_s, iso_limits[i]);
        if (rms_c > iso_limits[i]) printf("  ← EXCEEDS LIMIT");
        printf("\n");
    }

    /* Suspension travel analysis */
    printf("\n=== Suspension Travel (Rattle Space) ===\n");
    SecondOrderSystem travel_c, travel_s;
    quarter_car_suspension_travel(&comfort_tune, &travel_c);
    quarter_car_suspension_travel(&sport_tune, &travel_s);

    printf("Comfort tune travel DC gain: %.4f m/m (%.1f mm/cm bump)\n",
           travel_c.K, travel_c.K * 1000.0);
    printf("Sport tune travel DC gain:   %.4f m/m (%.1f mm/cm bump)\n",
           travel_s.K, travel_s.K * 1000.0);

    /* Design recommendation */
    printf("\n=== Design Summary ===\n");
    printf("Comfort tune:  Body ζ=%.3f (target 0.2-0.4) — %s\n",
           body_comfort.zeta,
           (body_comfort.zeta >= 0.2 && body_comfort.zeta <= 0.4)
           ? "✓ IN RANGE" : "✗ OUT OF RANGE");
    printf("Sport tune:    Body ζ=%.3f (target 0.2-0.4) — %s\n",
           body_sport.zeta,
           (body_sport.zeta >= 0.2 && body_sport.zeta <= 0.4)
           ? "✓ IN RANGE" : "✗ OUT OF RANGE");
    printf("Wheel hop damping should be > 0.2 to prevent wheel bounce.\n");
    printf("Comfort: wheel ζ=%.3f, Sport: wheel ζ=%.3f\n",
           wheel_comfort.zeta, wheel_sport.zeta);

    printf("\n=== Analysis Complete ===\n");
    return 0;
}
