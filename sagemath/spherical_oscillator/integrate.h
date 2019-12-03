long double get_sigma(charge q, coordinate r);
long double get_dS_dtheta(coordinate r, angle theta);
int calc_R_lw(charge q, timevalue t, coordinate R0, coordinate r0, velocity v0, acceleration a0, angle theta, timevalue * pt_zap, coordinate * pr_zap, distance * pR_zap, coordinate r_min, distance * R_lw_zap);
long double get_cos_alpha_zap(coordinate R0, angle theta, coordinate r_zap, distance R_zap);
long double get_aR_zap(coordinate R0, angle theta, coordinate r_zap, acceleration a_zap);
field get_E_minus_grad_phi_R0(angle theta, velocity v_zap, distance R_zap, long double aR_zap, distance R_lw_zap, long double cos_alpha_zap);
field get_E_minus_1_c_dA_dt_R0(angle theta, velocity v_zap, acceleration a_zap, distance R_zap, long double aR_zap, distance R_lw_zap);
int integral_phi(charge q, timevalue t, coordinate R0, coordinate r0, velocity v0, double a0, coordinate r_min, potential * result);
int integral_phi_and_E(charge q, timevalue t, coordinate R0, coordinate r0, velocity v0, acceleration a0, field * pE_minus_grad_phi_R0, field *pE_minus_1_c_dA_dt_R0, coordinate r_min, potential *phi);

