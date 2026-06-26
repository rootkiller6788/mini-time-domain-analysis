/* Example 3: PID Controller Design for Type 2 Performance
 * Plant Gp=2/(3s+1) Type 0, PI: Gc=Kp+Ki/s => Type 1
 * Desired: Kv=10 via PI tuning */
#include "steady_state_error.h"
#include "error_constants.h"
#include "system_type.h"
#include <stdio.h>
#include <stdlib.h>
int main(void){
    double plant_dc=2.0, desired_Kv=10.0, Kp, Ki;
    st_design_pi_for_type1(desired_Kv, plant_dc, &Kp, &Ki);
    printf("PI Controller Design for Kv=%.1f:
", desired_Kv);
    printf("  Plant DC gain = %.1f
", plant_dc);
    printf("  Kp = %.4f, Ki = %.4f
", Kp, Ki);
    printf("  Zero location: s = -Ki/Kp = %.4f
", -Ki/Kp);
    printf("
Open-loop: G_ol = (%.4f + %.4f/s) * 2/(3s+1)
", Kp, Ki);
    printf("Kv = lim s*G_ol = Ki*2 = %.4f (target: %.1f)
", Ki*2, desired_Kv);
    printf("e_ss_ramp = 1/Kv = %.4f
", 1.0/(Ki*2));
    printf("e_ss_step = 0 (Type 1+)
");
    return 0;
}