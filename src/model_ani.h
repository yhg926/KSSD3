#ifndef MODEL_ANI_H
#define MODEL_ANI_H

// model 1: 3-way linear model parameters (NUM_CODENS==9 only)

// optimized parametos for 3 denominator of 3-way linear model
double opt_denom_params4x3[12] =
    {
        1.0163855, 0.9824051, 1.0490899, 1.0110220, // denom1
        1.0196821, 0.9817771, 1.0463704, 1.0244149, // denom2
        1.0074051, 1.0018403, 0.9465413, 0.9405132  // denom3
};

//  3-way linear model coefficients

double linear_coeffs_3way_9CODENs[17] =
    {
        10.7394186713621,      // (Intercept)
        -2.92835942110571e-08, // XnY_ctx
        1.66112325921789e-05,  // N_diff_obj_section
        -3.48757817176062e-05, // N_mut2_ctx
        -1.82725607905918e-05, // N_diff_obj
        -1300.37798794877,     // denom1
        3151.11538582057,      // denom2
        -1824.9657489102,      // denom3
        25441.8492503678,      // XnY_ctx:denom1
        35763.8363137582,      // N_diff_obj_section:denom1
        -3901.01958179381,     // N_mut2_ctx:denom1
        -21820.9227249391,     // XnY_ctx:denom2
        -30943.0950958969,     // N_diff_obj_section:denom2
        10100.3465908889,      // N_mut2_ctx:denom2
        -3669.67487948577,     // XnY_ctx:denom3
        -4909.59256128814,     // N_diff_obj_section:denom3
        -5653.59698016275,     // N_mut2_ctx:denom3
};

typedef struct
{
    uint32_t XnY_ctx;
    uint32_t N_diff_obj_section;
    uint32_t N_mut2_ctx;
    uint32_t N_diff_obj;
} ani_features_t;

#define EPSILON (1e-8)
static inline double lm3ways_dist_from_features(ani_features_t *features)
{
    // compute 3 denomonators
    double denom1 = 1 / (opt_denom_params4x3[0] * features->XnY_ctx + opt_denom_params4x3[1] * features->N_diff_obj_section + opt_denom_params4x3[2] * features->N_mut2_ctx + opt_denom_params4x3[3] * features->N_diff_obj + EPSILON);
    double denom2 = 1 / (opt_denom_params4x3[4] * features->XnY_ctx + opt_denom_params4x3[5] * features->N_diff_obj_section + opt_denom_params4x3[6] * features->N_mut2_ctx + opt_denom_params4x3[7] * features->N_diff_obj + EPSILON);
    double denom3 = 1 / (opt_denom_params4x3[8] * features->XnY_ctx + opt_denom_params4x3[9] * features->N_diff_obj_section + opt_denom_params4x3[10] * features->N_mut2_ctx + opt_denom_params4x3[11] * features->N_diff_obj + EPSILON);

    // compute 3-way linear model
    double dist = linear_coeffs_3way_9CODENs[0] +
                  linear_coeffs_3way_9CODENs[1] * features->XnY_ctx +
                  linear_coeffs_3way_9CODENs[2] * features->N_diff_obj_section +
                  linear_coeffs_3way_9CODENs[3] * features->N_mut2_ctx +
                  linear_coeffs_3way_9CODENs[4] * features->N_diff_obj +
                  linear_coeffs_3way_9CODENs[5] * denom1 +
                  linear_coeffs_3way_9CODENs[6] * denom2 +
                  linear_coeffs_3way_9CODENs[7] * denom3 +
                  linear_coeffs_3way_9CODENs[8] * features->XnY_ctx * denom1 +
                  linear_coeffs_3way_9CODENs[9] * features->N_diff_obj_section * denom1 +
                  linear_coeffs_3way_9CODENs[10] * features->N_mut2_ctx * denom1 +
                  linear_coeffs_3way_9CODENs[11] * features->XnY_ctx * denom2 +
                  linear_coeffs_3way_9CODENs[12] * features->N_diff_obj_section * denom2 +
                  linear_coeffs_3way_9CODENs[13] * features->N_mut2_ctx * denom2 +
                  linear_coeffs_3way_9CODENs[14] * features->XnY_ctx * denom3 +
                  linear_coeffs_3way_9CODENs[15] * features->N_diff_obj_section * denom3 +
                  linear_coeffs_3way_9CODENs[16] * features->N_mut2_ctx * denom3;
    if(dist < 0) dist = 0; 
    return dist;
}
#endif