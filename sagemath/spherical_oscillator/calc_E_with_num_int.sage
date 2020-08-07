attach("../field_of_deyna_cylinder/num_int.sage")

def phi_and_E_integrand_phi(q, t, R0, r0, v0, a0, theta, r_min):
    (phi, A, E1, E2, E, error) = phi_and_E_integrand(q, t, R0, r0, v0, a0, theta, r_min)
    return phi
def phi_and_E_integrand_A(q, t, R0, r0, v0, a0, theta, r_min):
    (phi, A, E1, E2, E, error) = phi_and_E_integrand(q, t, R0, r0, v0, a0, theta, r_min)
    return A

def phi_and_E_integrand_E(q, t, R0, r0, v0, a0, theta, r_min):
    (phi, A, E1, E2, E, error) = phi_and_E_integrand(q, t, R0, r0, v0, a0, theta, r_min)
    return E
def phi_and_E_integrand_E1(q, t, R0, r0, v0, a0, theta, r_min):
    (phi, A, E1, E2, E, error) = phi_and_E_integrand(q, t, R0, r0, v0, a0, theta, r_min)
    return E1
def phi_and_E_integrand_E2(q, t, R0, r0, v0, a0, theta, r_min):
    (phi, A, E1, E2, E, error) = phi_and_E_integrand(q, t, R0, r0, v0, a0, theta, r_min)
    return E2

calc_phi = lambda q, t, R0, r0, v0, a0, theta, r_min : phi_and_E_integrand_phi (q, t, R0, r0, v0, a0, theta, r_min)
calc_A   = lambda q, t, R0, r0, v0, a0, theta, r_min : phi_and_E_integrand_A   (q, t, R0, r0, v0, a0, theta, r_min)
calc_E   = lambda q, t, R0, r0, v0, a0, theta, r_min : phi_and_E_integrand_E   (q, t, R0, r0, v0, a0, theta, r_min)
calc_E1  = lambda q, t, R0, r0, v0, a0, theta, r_min : phi_and_E_integrand_E1  (q, t, R0, r0, v0, a0, theta, r_min)
calc_E2  = lambda q, t, R0, r0, v0, a0, theta, r_min : phi_and_E_integrand_E2  (q, t, R0, r0, v0, a0, theta, r_min)

k_E     = lambda q, t, R0, r0, v0, a0, r_min : num_int(lambda theta : calc_E   ( q, t, R0, r0, v0, a0, theta, r_min) / m * sin(theta) * R0^2, 0, pi )
k_E1    = lambda q, t, R0, r0, v0, a0, r_min : num_int(lambda theta : calc_E1  ( q, t, R0, r0, v0, a0, theta, r_min) / m * sin(theta) * R0^2, 0, pi )
k_E2    = lambda q, t, R0, r0, v0, a0, r_min : num_int(lambda theta : calc_E2  ( q, t, R0, r0, v0, a0, theta, r_min) / m * sin(theta) * R0^2, 0, pi )
int_phi = lambda q, t, R0, r0, v0, a0, r_min : num_int(lambda theta : calc_phi ( q, t, R0, r0, v0, a0, theta, r_min) * 2 * pi * sin(theta) * r0^2, 0, pi ) * q / (4 * pi.n() * r0^2)
int_A   = lambda q, t, R0, r0, v0, a0, r_min : num_int(lambda theta : calc_A   ( q, t, R0, r0, v0, a0, theta, r_min) * 2 * pi * sin(theta) * r0^2, 0, pi ) * q / (4 * pi.n() * r0^2)
int_E   = lambda q, t, R0, r0, v0, a0, r_min : num_int(lambda theta : calc_E   ( q, t, R0, r0, v0, a0, theta, r_min) * 2 * pi * sin(theta) * r0^2, 0, pi ) * q / (4 * pi.n() * r0^2)
int_E1  = lambda q, t, R0, r0, v0, a0, r_min : num_int(lambda theta : calc_E1  ( q, t, R0, r0, v0, a0, theta, r_min) * 2 * pi * sin(theta) * r0^2, 0, pi ) * q / (4 * pi.n() * r0^2)
int_E2  = lambda q, t, R0, r0, v0, a0, r_min : num_int(lambda theta : calc_E2  ( q, t, R0, r0, v0, a0, theta, r_min) * 2 * pi * sin(theta) * r0^2, 0, pi ) * q / (4 * pi.n() * r0^2)
