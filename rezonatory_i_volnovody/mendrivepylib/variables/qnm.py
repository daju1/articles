from sage.calculus.calculus import var
from sage.symbolic.assumptions import assume

L_z =  var("L_z", domain='real')        # длина резонатора вдоль z, см
omega = var('omega', domain='complex')  # частота в общем виде
                                        # применяется до подстановки
                                        # для решения в общем виде обеих задач: QNM и drive
omega_drive = var('omega_drive', domain='real') # частота возбуждения
omega_re = var("omega_re", domain='real')       # компоненты комплексной частоты в задаче QNM
omega_im = var("omega_im", domain='real')       # компоненты комплексной частоты в задаче QNM

# ======================================================
# 1. Продольное волновое число k_z (вещественное при резонансе)
# ======================================================
k_z = var('k_z', domain='real')

k_y = var ('k_y')
assume(k_y, 'complex')
