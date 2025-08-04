#ifndef MODEL_ANI_H
#define MODEL_ANI_H

// model 1: 3-way linear model parameters (NUM_CODENS==9 only)

// optimized parameters for 3 denominator of 3-way linear model
static const double opt_denom_params4x3[12] =
    {
        1.0090536851341, 1.0013604061893, 1.03540226374509, 1.04270754727088, // denom1
        1.0108940845453, 1.0010797461052, 0.98914550652945, 1.01536320288592, // denom2
        1.0061957741437, 0.9957814716168, 0.98352379677340, 1.00951820561544  // denom3
};
// ANI>95 subset parameters
static const double N95opt_denom_params4x3[12] =
    {
        1.0196821174688, 0.9817771341602, 1.04637038371202, 1.02441485561367, // denom1
        1.0163855061091, 0.9824051281495, 1.04908985213703, 1.01102197270052, // denom2
        1.0074050784795, 1.0018402553389, 0.94654127396428, 0.94051319172603  // denom3
};

//  3-way linear model coefficients
static const double linear_coeffs_3way_9CODENs[17] =
    {
        -0.587100397752123,    // (Intercept)
        -7.93027271403258e-08, // XnY_ctx
        4.98481429879899e-06,  // N_diff_obj_section
        -1.63903694463017e-05, // N_mut2_ctx
        -4.93356600478016e-06, // N_diff_obj
        -14843.5350131617,     // denom1
        306478.144745683,      // denom2
        -290254.345158429,     // denom3
        -1781.00695309133,     // XnY_ctx:denom1
        -1091.99589131905,     // N_diff_obj_section:denom1
        327.6575530681,        // N_mut2_ctx:denom1
        47395.255680747,       // XnY_ctx:denom2
        47798.7112013701,      // N_diff_obj_section:denom2
        -35574.6303634354,     // N_mut2_ctx:denom2
        -45398.423188085,      // XnY_ctx:denom3
        -46461.4920783106,     // N_diff_obj_section:denom3
        35045.3583333325,      // N_mut2_ctx:denom3
};
// ANI>95 subset of 3-way linear model coefficients
static const double N95linear_coeffs_3way_9CODENs[17] =
    {
        10.7394186800999,      // (Intercept)
        -2.92835943410177e-08, // XnY_ctx
        1.66112325940088e-05,  // N_diff_obj_section
        -3.48757817098253e-05, // N_mut2_ctx
        -1.82725607919882e-05, // N_diff_obj
        3151.11534599359,      // denom1
        -1300.37794535582,     // denom2
        -1824.96575177919,     // denom3
        -21820.9227396359,     // XnY_ctx:denom1
        -30943.09511694,       // N_diff_obj_section:denom1
        10100.3465983902,      // N_mut2_ctx:denom1
        25441.8492675468,      // XnY_ctx:denom2
        35763.8363380932,      // N_diff_obj_section:denom2
        -3901.01958522226,     // N_mut2_ctx:denom2
        -3669.6748819999,      // XnY_ctx:denom3
        -4909.59256464219,     // N_diff_obj_section:denom3
        -5653.59698386914,     // N_mut2_ctx:denom3

};

typedef struct
{
    uint32_t XnY_ctx;
    uint32_t N_diff_obj_section;
    uint32_t N_mut2_ctx;
    uint32_t N_diff_obj;
} ani_features_t;

#define EPSILON (1e-8)
static inline double lm3ways_dist_from_features_core(ani_features_t *features, const double p[12], const double coeffs[17])
{
    // compute 3 denomonators
    double denom1 = 1 / (p[0] * features->XnY_ctx + p[1] * features->N_diff_obj_section + p[2] * features->N_mut2_ctx + p[3] * features->N_diff_obj + EPSILON);
    double denom2 = 1 / (p[4] * features->XnY_ctx + p[5] * features->N_diff_obj_section + p[6] * features->N_mut2_ctx + p[7] * features->N_diff_obj + EPSILON);
    double denom3 = 1 / (p[8] * features->XnY_ctx + p[9] * features->N_diff_obj_section + p[10] * features->N_mut2_ctx + p[11] * features->N_diff_obj + EPSILON);
    // compute 3-way linear model
    double dist = coeffs[0] +
                  coeffs[1] * features->XnY_ctx +
                  coeffs[2] * features->N_diff_obj_section +
                  coeffs[3] * features->N_mut2_ctx +
                  coeffs[4] * features->N_diff_obj +
                  coeffs[5] * denom1 +
                  coeffs[6] * denom2 +
                  coeffs[7] * denom3 +
                  coeffs[8] * features->XnY_ctx * denom1 +
                  coeffs[9] * features->N_diff_obj_section * denom1 +
                  coeffs[10] * features->N_mut2_ctx * denom1 +
                  coeffs[11] * features->XnY_ctx * denom2 +
                  coeffs[12] * features->N_diff_obj_section * denom2 +
                  coeffs[13] * features->N_mut2_ctx * denom2 +
                  coeffs[14] * features->XnY_ctx * denom3 +
                  coeffs[15] * features->N_diff_obj_section * denom3 +
                  coeffs[16] * features->N_mut2_ctx * denom3;
    return dist;
}

// 2. 3-ways linear model distance
static inline double lm3ways_dist_from_features(ani_features_t *features)
{
    if (features-> XnY_ctx == 0) return 1; // if no ctx, return 1
    else if (features->N_diff_obj == 0) return 0; // if no diff obj, return 0
    // use optimized parameters lm prediction     
    double dist = lm3ways_dist_from_features_core(features, opt_denom_params4x3, linear_coeffs_3way_9CODENs);
    if (dist < 0.05) dist = lm3ways_dist_from_features_core(features, N95opt_denom_params4x3, N95linear_coeffs_3way_9CODENs);
    if (dist < 0) dist = 0;
    return dist;
}


#endif