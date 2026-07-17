from sage.calculus.calculus import var
from sage.symbolic.assumptions import assume

mu_l_zy_c = var('mu_l_zy_c')
mu_l_yz_c = var('mu_l_yz_c')

mu_r_xx_c = var('mu_r_xx_c')
mu_r_yy_c = var('mu_r_yy_c')
mu_r_zz_c = var('mu_r_zz_c')

mu_l_xx_c = var('mu_l_xx_c')
mu_l_yy_c = var('mu_l_yy_c')
mu_l_zz_c = var('mu_l_zz_c')

assume(mu_l_xx_c, 'complex')
assume(mu_l_yy_c, 'complex')
assume(mu_l_zz_c, 'complex')

mu_l_xx = var('mu_l_xx')
mu_l_yy = var('mu_l_yy')
mu_l_zz = var('mu_l_zz')

assume(mu_l_xx, 'real')
assume(mu_l_yy, 'real')
assume(mu_l_zz, 'real')

assume(mu_l_xx>0)
assume(mu_l_yy>0)
assume(mu_l_zz>0)

mu_r_xx = var('mu_r_xx')
mu_r_yy = var('mu_r_yy')
mu_r_zz = var('mu_r_zz')

assume(mu_r_xx, 'real')
assume(mu_r_yy, 'real')
assume(mu_r_zz, 'real')

assume(mu_r_xx>0)
assume(mu_r_yy>0)
assume(mu_r_zz>0)

mu_r_zy_c = var('mu_r_zy_c')
mu_r_yz_c = var('mu_r_yz_c')

epsilon_l_xx_c = var('epsilon_l_xx_c')
epsilon_l_yy_c = var('epsilon_l_yy_c')
epsilon_l_zz_c = var('epsilon_l_zz_c')

assume(epsilon_l_xx_c, 'complex')
assume(epsilon_l_yy_c, 'complex')
assume(epsilon_l_zz_c, 'complex')

epsilon_l_xx = var('epsilon_l_xx')
epsilon_l_yy = var('epsilon_l_yy')
epsilon_l_zz = var('epsilon_l_zz')

assume(epsilon_l_xx, 'real')
assume(epsilon_l_yy, 'real')
assume(epsilon_l_zz, 'real')

assume(epsilon_l_xx>0)
assume(epsilon_l_yy>0)
assume(epsilon_l_zz>0)

epsilon_r_xx_c = var('epsilon_r_xx_c')
epsilon_r_yy_c = var('epsilon_r_yy_c')
epsilon_r_zz_c = var('epsilon_r_zz_c')

assume(epsilon_r_xx_c, 'complex')
assume(epsilon_r_yy_c, 'complex')
assume(epsilon_r_zz_c, 'complex')

epsilon_r_xx = var('epsilon_r_xx')
epsilon_r_yy = var('epsilon_r_yy')
epsilon_r_zz = var('epsilon_r_zz')

assume(epsilon_r_xx, 'real')
assume(epsilon_r_yy, 'real')
assume(epsilon_r_zz, 'real')

assume(epsilon_r_xx>0)
assume(epsilon_r_yy>0)
assume(epsilon_r_zz>0)

epsilon_r_xy_c = var('epsilon_r_xy_c')
epsilon_r_yx_c = var('epsilon_r_yx_c')

epsilon_r_xz_c = var('epsilon_r_xz_c')
epsilon_r_zx_c = var('epsilon_r_zx_c')

epsilon_r_zy_c = var('epsilon_r_zy_c')
epsilon_r_yz_c = var('epsilon_r_yz_c')

mu_0 = var("mu_0")
assume(mu_0, "real")
assume(mu_0>0)

epsilon_0 = var("epsilon_0")
assume(epsilon_0, "real")
assume(epsilon_0>0)

# sigma_0 = var("sigma_0")
# assume(sigma_0, "real")
# assume(sigma_0>0)

# left conductor
sigma_e_l_xx = var('sigma_e_l_xx')
sigma_e_l_yy = var('sigma_e_l_yy')
sigma_e_l_zz = var('sigma_e_l_zz')

assume(sigma_e_l_xx, 'real')
assume(sigma_e_l_yy, 'real')
assume(sigma_e_l_zz, 'real')

assume(sigma_e_l_xx>0)
assume(sigma_e_l_yy>0)
assume(sigma_e_l_zz>0)

# left conductor
sigma_m_l_xx = var('sigma_m_l_xx')
sigma_m_l_yy = var('sigma_m_l_yy')
sigma_m_l_zz = var('sigma_m_l_zz')

assume(sigma_m_l_xx, 'real')
assume(sigma_m_l_yy, 'real')
assume(sigma_m_l_zz, 'real')

assume(sigma_m_l_xx>0)
assume(sigma_m_l_yy>0)
assume(sigma_m_l_zz>0)

#right conductor
sigma_e_r_xx = var('sigma_e_r_xx')
sigma_e_r_yy = var('sigma_e_r_yy')
sigma_e_r_zz = var('sigma_e_r_zz')

assume(sigma_e_r_xx, 'real')
assume(sigma_e_r_yy, 'real')
assume(sigma_e_r_zz, 'real')

assume(sigma_e_r_xx>0)
assume(sigma_e_r_yy>0)
assume(sigma_e_r_zz>0)

#right conductor
sigma_m_r_xx = var('sigma_m_r_xx')
sigma_m_r_yy = var('sigma_m_r_yy')
sigma_m_r_zz = var('sigma_m_r_zz')

assume(sigma_m_r_xx, 'real')
assume(sigma_m_r_yy, 'real')
assume(sigma_m_r_zz, 'real')

assume(sigma_m_r_xx>0)
assume(sigma_m_r_yy>0)
assume(sigma_m_r_zz>0)

kappa_vacuum = var ('kappa_vacuum')
assume(kappa_vacuum, 'real')
assume(kappa_vacuum>0)

# left conductor
kappa_l = var('kappa_l')
assume(kappa_l, 'complex')

# right conductor
kappa_r = var('kappa_r')
assume(kappa_r, 'complex')

k_y_E = var ('k_y_E')
assume(k_y_E, 'real')

k_y_H = var ('k_y_H')
assume(k_y_H, 'real')

a = var('a')
assume(a, 'real')
assume(a>0)

x = var('x')
y = var('y')
z = var('z')

assume(x, 'real')
assume(y, 'real')
assume(z, 'real')

b = var('b')
m = var('m')

c = var('c')
assume(c, 'real')
assume(c>0)
