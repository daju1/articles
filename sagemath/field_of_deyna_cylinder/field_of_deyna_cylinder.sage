import sys
reload(sys)
sys.setdefaultencoding('utf8')
import numpy as np


attach("float_formatting.sage")


def get_integrand_view(f):
    return f(x)


def my_numerical_integral(f, a, b):
    from scipy import integrate

    log_fn = 'field_of_deyna_cylinder_my_numerical_integral.txt'
    # log_fn = 'field_of_deyna_cylinder.txt'

    # import traceback
    # traceback.print_stack()

    to_call_integration = True

    if None == f:
        to_call_integration = False

    import inspect
    stack = inspect.stack()
    for frame in stack:
        func_name = frame[3]
        print "func_name = ", func_name
        if ('get_integrand_view' == func_name):
            to_call_integration = False
            break;

    print "stack = ", stack
    print "f = ", f
    print "integrand = ", get_integrand_view(f)
    print "a = ", a
    print "b = ", b
    print "to_call_integration = ", to_call_integration

    # file = open('field_of_deyna_cylinder_my_numerical_integral.txt', 'a')
    file = open(log_fn, 'a')
    file.write('\n')
    file.write("stack = " + str(stack))
    file.write('\n')
    file.write("f = " + str(f))
    file.write('\n')
    file.write("integrand = " + str(get_integrand_view(f)))
    file.write('\n')
    file.write("a = " + str(a) + ", b = " + str(b))
    file.write('\n')
    file.write("to_call_integration = " + str(to_call_integration))
    file.write('\n\n')
    file.close()

    if not to_call_integration:
        return None

    try:
        integral = integrate.quad(f, a, b)
        # integral = numerical_integral(f, a, b)

        print "integral = ", integral

        # file = open(log_fn, 'a')
        file = open(log_fn, 'a')
        file.write('\n')
        file.write("integral = " + str(integral))
        file.write('\n\n')
        file.close()

        result = integral[0]
        return result

    except Exception as ex:

        print "Exception ex = ", str(ex)
        print "f = ", f
        print "integrand = ", get_integrand_view(f)
        print "a = ", a
        print "b = ", b

        file = open(log_fn, 'a')
        file.write('\n')
        file.write("Exception ex = " + str(ex))
        file.write('\n')
        file.write("f = " + str(f))
        file.write('\n')
        file.write("integrand = " + str(get_integrand_view(f)))
        file.write('\n')
        file.write("a = " + str(a) + ", b = " + str(b))
        file.write('\n\n')
        file.write('\n\n')
        file.close()

        if 'unable to simplify to float approximation' == str(ex):
            raise ex

        integral = numerical_integral(f, a, b)

        print "integral = ", integral

        file = open(log_fn, 'a')
        file.write('\n')
        file.write("integral = " + str(integral))
        file.write('\n\n')
        file.close()

        result = integral[0]
        print "result = ", result
        return result


zj = var("zj")
za = var("za")
rj = var("rj")
ra = var("ra")
phi = var("phi")


assume(rj>0)
assume(ra>0)

# integral(1/sqrt(rj^2+ra^2+(zj-za)^2-2*rj*ra*cos(phi)), (phi, 0, 2*pi))
# Ivarphi__j = integral(1/sqrt(rho__j^2+rho__a^2+(z__j-z__a)^2-2*rho__j*rho__a*cos(varphi__j)), (varphi__j,0,2*pi))
# RuntimeError: ECL says: Error executing code in Maxima:

# I_phi_j = integral(1/sqrt(rj^2+ra^2+(zj-za)^2-2*rj*ra*cos(phi)), phi)
# print "I_phi_j =", I_phi_j
# SageMath version 8.4, Release Date: 2018-10-17
# BUG
# I_phi_j = 1/2*log(cos(phi)^2 + sin(phi)^2 + 2*sin(phi) + 1) - 1/2*log(cos(phi)^2 + sin(phi)^2 - 2*sin(phi) + 1)

rja2 = (rj-ra)^2+(zj-za)^2
module = - 4*rj*ra / rja2

I_phi_j = 2 * elliptic_f(phi/2, module) / sqrt(rja2)
print "I_phi_j =", I_phi_j
I_phi_j_diff_phi = I_phi_j.diff(phi)
print "I_phi_j_diff_phi =", I_phi_j_diff_phi
# I_phi_j.diff(phi) = 1/(sqrt(4*ra*rj*sin(1/2*phi)^2/((ra - rj)^2 + (za - zj)^2) + 1)*sqrt((ra - rj)^2 + (za - zj)^2))
# sage: I_phi_j.diff(phi)
#                                 1                                
# -----------------------------------------------------------------
#        _____________________________                             
#       /               2/phi\                                     
#      /     4*ra*rj*sin |---|            _________________________
#     /                  \ 2 /           /          2            2 
#    /    ----------------------- + 1 *\/  (ra - rj)  + (za - zj)  
#   /              2            2                                  
# \/      (ra - rj)  + (za - zj)                                   

# I_phi_j.diff(phi) = 1/(sqrt(4*ra*rj*sin(1/2*phi)^2/((ra - rj)^2 + (za - zj)^2) + 1)*sqrt((ra - rj)^2 + (za - zj)^2))

I_phi_j_diff_phi_substituted = I_phi_j_diff_phi.substitute((sin(1/2*phi)^2)==(1-cos(phi))/2)
print "I_phi_j_diff_phi_substituted =", I_phi_j_diff_phi_substituted
# I_phi_j.diff(phi) = 1/(sqrt(4*ra*rj*sin(1/2*phi)^2/((ra - rj)^2 + (za - zj)^2) + 1)*sqrt((ra - rj)^2 + (za - zj)^2))
# sage: (I_phi_j.diff(phi).substitute((sin(1/2*phi)^2)==(1-cos(phi))/2))
#                                 1                                
# -----------------------------------------------------------------
#      _______________________________    _________________________
#     /    2*ra*rj*(cos(phi) - 1)        /          2            2 
#    /  - ----------------------- + 1 *\/  (ra - rj)  + (za - zj)  
#   /              2            2                                  
# \/      (ra - rj)  + (za - zj)                                   

print "I_phi_j_diff_phi_substituted.full_simplify() =", I_phi_j_diff_phi_substituted.full_simplify()
# I_phi_j_diff_phi_substituted.full_simplify() = 1/(sqrt(ra^2 - 2*ra*rj + rj^2 + za^2 - 2*za*zj + zj^2)*sqrt(-(2*ra*rj*cos(phi) - ra^2 - rj^2 - za^2 + 2*za*zj - zj^2)/(ra^2 - 2*ra*rj + rj^2 + za^2 - 2*za*zj + zj^2)))

print "sqrt(I_phi_j_diff_phi_substituted^2) =", sqrt(I_phi_j_diff_phi_substituted^2)
# sqrt(I_phi_j_diff_phi_substituted^2) = sqrt(-1/(((ra - rj)^2 + (za - zj)^2)*(2*ra*rj*(cos(phi) - 1)/((ra - rj)^2 + (za - zj)^2) - 1)))

print "sqrt(I_phi_j_diff_phi_substituted^2).full_simplify() =", sqrt(I_phi_j_diff_phi_substituted^2).full_simplify()
# sqrt(I_phi_j_diff_phi_substituted^2).full_simplify() = sqrt(-1/(2*ra*rj*cos(phi) - ra^2 - rj^2 - za^2 + 2*za*zj - zj^2))
# sage: sqrt(I_phi_j_diff_phi_substituted^2).full_simplify()
#     ______________________________________________________
#    /                         -1                           
#   /  ---------------------------------------------------- 
#  /       2                        2     2               2 
#\/    - ra  + 2*ra*rj*cos(phi) - rj  - za  + 2*za*zj - zj  


print "elliptic_kc(0.5)=", elliptic_kc(0.5)
# 1.85407467730137

print "elliptic_kc(-0.5)=", elliptic_kc(-0.5)
# 1.41573720842596

m = var("m")
print elliptic_kc(m).diff(m)
# -1/2*((m - 1)*elliptic_kc(m) + elliptic_ec(m))/((m - 1)*m)


Iphi=4*elliptic_kc(module) / sqrt(rja2)
print "Iphi =", Iphi
# Iphi = 4*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (za - zj)^2))/sqrt((ra - rj)^2 + (za - zj)^2)
# sage: %display ascii_art
# sage: Iphi
#    /       -4*ra*rj        \
# 4*K|-----------------------|
#    |         2            2|
#    \(ra - rj)  + (za - zj) /
# ----------------------------
#    _________________________
#   /          2            2 
# \/  (ra - rj)  + (za - zj)  

# from sage.symbolic.integration.integral import definite_integral
# II_phi = definite_integral(1/sqrt(rj^2+ra^2+(zj-za)^2-2*rj*ra*cos(phi)), phi,0,2*pi)
# print "II_phi =", II_phi

print "Iphi.substitute(rj=1, zj=1, ra=2, za=0) =", Iphi.substitute(rj=1, zj=1, ra=2, za=0).n()
# 2.85516399176746

print "Iphi.substitute(rj=1, zj=1, ra=1, za=-1)", Iphi.substitute(rj=1, zj=1, ra=1, za=-1).n()
# 2.62205755429212

dIphi_dza = Iphi.diff(za)
print "dIphi_dza =", dIphi_dza
# dIphi_dza = -4*(za - zj)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (za - zj)^2))/((ra - rj)^2 + (za - zj)^2)^(3/2) + 4*((4*ra*rj/((ra - rj)^2 + (za - zj)^2) + 1)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (za - zj)^2)) - elliptic_ec(-4*ra*rj/((ra - rj)^2 + (za - zj)^2)))*(za - zj)/(((ra - rj)^2 + (za - zj)^2)^(3/2)*(4*ra*rj/((ra - rj)^2 + (za - zj)^2) + 1))

#             //        4*ra*rj            \  /       -4*ra*rj        \    /       -4*ra*rj        \\                /       -4*ra*rj        \
# 4*(za - zj)*||----------------------- + 1|*K|-----------------------| - E|-----------------------||   4*(za - zj)*K|-----------------------|
#             ||         2            2    |  |         2            2|    |         2            2||                |         2            2|
#             \\(ra - rj)  + (za - zj)     /  \(ra - rj)  + (za - zj) /    \(ra - rj)  + (za - zj) //                \(ra - rj)  + (za - zj) /
# --------------------------------------------------------------------------------------------------- - --------------------------------------
#                                                                             3/2                                                     3/2     
#                      /        4*ra*rj            \ /         2            2\                               /         2            2\        
#                      |----------------------- + 1|*\(ra - rj)  + (za - zj) /                               \(ra - rj)  + (za - zj) /        
#                      |         2            2    |                            
#                      \(ra - rj)  + (za - zj)     /                            

print "dIphi_dza.full_simplify() =", dIphi_dza.full_simplify()


dIphi_dra = Iphi.diff(ra)
print "dIphi_dra =", dIphi_dra
# dIphi_dra = -4*(ra - rj)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (za - zj)^2))/((ra - rj)^2 + (za - zj)^2)^(3/2) + 2*sqrt((ra - rj)^2 + (za - zj)^2)*((4*ra*rj/((ra - rj)^2 + (za - zj)^2) + 1)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (za - zj)^2)) - elliptic_ec(-4*ra*rj/((ra - rj)^2 + (za - zj)^2)))*(2*(ra - rj)*ra*rj/((ra - rj)^2 + (za - zj)^2)^2 - rj/((ra - rj)^2 + (za - zj)^2))/(ra*(4*ra*rj/((ra - rj)^2 + (za - zj)^2) + 1)*rj)

                                                                                
#                                              //        4*ra*rj            \  /       -4*ra*rj        \    /       -4*ra*rj        \\ /    2*ra*rj*(ra - rj) 
#                /       -4*ra*rj        \   2*||----------------------- + 1|*K|-----------------------| - E|-----------------------||*|----------------------
#   4*(ra - rj)*K|-----------------------|     ||         2            2    |  |         2            2|    |         2            2|| |                      
#                |         2            2|     \\(ra - rj)  + (za - zj)     /  \(ra - rj)  + (za - zj) /    \(ra - rj)  + (za - zj) // |/         2           
#                \(ra - rj)  + (za - zj) /                                                                                             \\(ra - rj)  + (za - zj
# - -------------------------------------- + -----------------------------------------------------------------------------------------------------------------
#                                 3/2                                                                                   /        4*ra*rj            \        
#        /         2            2\                                                                                ra*rj*|----------------------- + 1|      
#        \(ra - rj)  + (za - zj) /                                                                                      |         2            2    |       
#                                                                                                                       \(ra - rj)  + (za - zj)     /         
#
#                                                                               
#                                    _________________________
# 
# ---- - -----------------------|*\/  (ra - rj)  + (za - zj)  
#    2            2            2|                             
#  2\    (ra - rj)  + (za - zj) |                             
# ) /                           /                             
# ------------------------------------------------------------
#                                                            
#                                                            
#                                                            
#       
#
if False:
    for r_j in np.arange(0.5, 5, 0.5):
        for r_a in np.arange(0.5, 5, 0.5):
            for z_j in np.arange(0.5, 5, 0.5):
                for z_a in np.arange(0.5, 5, 0.5):
                    if not ((r_j == r_a) and (z_j == z_a)):
                        print "dIphi_dza.substitute(rj=", r_j, ", zj=", z_j, ", ra=", r_a, ", za=", z_a, ") = ", dIphi_dza.substitute(rj=r_j, zj=z_j, ra=r_a, za=z_a).n()
                        #print "dIphi_dra.substitute(rj=", r_j, ", zj=", z_j, ", ra=", r_a, ", za=", z_a, ") = ", dIphi_dra.substitute(rj=r_j, zj=z_j, ra=r_a, za=z_a).n()

print "dIphi_dza.substitute(rj=1, zj=1, ra=1, za=-1) = ", dIphi_dza.substitute(rj=1, zj=1, ra=1, za=-1).n()
print "dIphi_dra.substitute(rj=1, zj=1, ra=1, za=-1) = ", dIphi_dra.substitute(rj=1, zj=1, ra=1, za=-1).n()
# dIphi_dza.substitute(rj=1, zj=1, ra=1, za=-1) =  0.955049447256928
# dIphi_dra.substitute(rj=1, zj=1, ra=1, za=-1) =  -0.355979329889132

# строго говоря размерность поверхностного тока должна быть приведена в соответствие с формулой js = c*[I x n]
# а размерность объёмного - формуле jv = c * rot(I)

cI0 = var("cI0")

full_volume_cylinder = False
full_volume_cylinder = True

if full_volume_cylinder:
    js = cI0
    jt = cI0
    jv = cI0 / rj
else:
    js = cI0 / rj
    jt = cI0 / ra
    jv = 0



zj1 = var("zj1")
zj2 = var("zj2")

rj1 = var("rj1")
rj2 = var("rj2")


Iphi_js_rj=Iphi*js*rj
print "Iphi_js_rj =", Iphi_js_rj

Iphi_jv_rj=Iphi*jv*rj
print "Iphi_jv_rj =", Iphi_jv_rj


Iphi_js_rj_diff_za = lambda c_I0, r_j, r_a, z_j, z_a : (Iphi_js_rj.diff(za)).substitute(ra==r_a).substitute(za==z_a).substitute(zj==z_j).substitute(rj==r_j).substitute(cI0==c_I0)
print "Iphi_js_rj_diff_za(cI0, rj, ra, zj, za) =", Iphi_js_rj_diff_za(cI0, rj, ra, zj, za)
print ""

file = open('field_of_deyna_cylinder.txt', 'a')
file.write("Iphi_js_rj_diff_za(cI0, rj, ra, zj, za) = " + str(Iphi_js_rj_diff_za(cI0, rj, ra, zj, za)))
file.write('\n\n')
file.close()

# At_diff_za = Iphi_js_rj_diff_za .substitute(zj==zj1) - Iphi_js_rj_diff_za .substitute(zj==zj2)
# At_diff_za = Iphi_js_rj_diff_za(cI0, rj, ra, zj, za).substitute(zj==zj1) - Iphi_js_rj_diff_za(cI0, rj, ra, zj, za).substitute(zj==zj2)
At_diff_za = lambda cI0, rj, ra, zj1, zj2, za : Iphi_js_rj_diff_za(cI0, rj, ra, zj1, za) - Iphi_js_rj_diff_za(cI0, rj, ra, zj2, za)
print "At_diff_za(cI0, rj, ra, zj1, zj2, za) =", At_diff_za(cI0, rj, ra, zj1, zj2, za)
print ""
# At_diff_za = 4*cI0*(z1 - za)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2))/((ra - rj)^2 + (z1 - za)^2)^(3/2) - 4*cI0*(z2 - za)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2))/((ra - rj)^2 + (z2 - za)^2)^(3/2) - 4*((4*ra*rj/((ra - rj)^2 + (z1 - za)^2) + 1)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2)) - elliptic_ec(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2)))*cI0*(z1 - za)/(((ra - rj)^2 + (z1 - za)^2)^(3/2)*(4*ra*rj/((ra - rj)^2 + (z1 - za)^2) + 1)) + 4*((4*ra*rj/((ra - rj)^2 + (z2 - za)^2) + 1)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2)) - elliptic_ec(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2)))*cI0*(z2 - za)/(((ra - rj)^2 + (z2 - za)^2)^(3/2)*(4*ra*rj/((ra - rj)^2 + (z2 - za)^2) + 1))

file = open('field_of_deyna_cylinder.txt', 'a')
file.write("At_diff_za(cI0, rj, ra, zj1, zj2, za) = " + str(At_diff_za(cI0, rj, ra, zj1, zj2, za)))
file.write('\n\n')
file.close()

Iphi_js_rj_diff_ra = lambda c_I0, r_j, r_a, z_j, z_a : (Iphi_js_rj.diff(ra)).substitute(ra==r_a).substitute(za==z_a).substitute(zj==z_j).substitute(rj==r_j).substitute(cI0==c_I0)
print "Iphi_js_rj_diff_ra(cI0, rj, ra, zj, za) =", Iphi_js_rj_diff_ra(cI0, rj, ra, zj, za)
print ""

file = open('field_of_deyna_cylinder.txt', 'a')
file.write("Iphi_js_rj_diff_ra(cI0, rj, ra, zj, za) = " + str(Iphi_js_rj_diff_ra(cI0, rj, ra, zj, za)))
file.write('\n\n')
file.close()

# As_diff_ra  =  - Iphi_js_rj_diff_ra.substitute(rj==rj1) + Iphi_js_rj_diff_ra.substitute(rj==rj2)
# As_diff_ra  =  - Iphi_js_rj_diff_ra(cI0, rj, ra, zj, za).substitute(rj==rj1) + Iphi_js_rj_diff_ra(cI0, rj, ra, zj, za).substitute(rj==rj2)
As_diff_ra  = lambda cI0, rj1, rj2, ra, zj, za : - Iphi_js_rj_diff_ra(cI0, rj1, ra, zj, za) + Iphi_js_rj_diff_ra(cI0, rj2, ra, zj, za)
print "As_diff_ra(cI0, rj1, rj2, ra, zj, za) =", As_diff_ra(cI0, rj1, rj2, ra, zj, za)
print ""
# As_diff_ra  = -4*cI0*(r1 - ra)*elliptic_kc(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2))/((r1 - ra)^2 + (za - zj)^2)^(3/2) + 4*cI0*(r2 - ra)*elliptic_kc(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2))/((r2 - ra)^2 + (za - zj)^2)^(3/2) + 2*sqrt((r1 - ra)^2 + (za - zj)^2)*((4*r1*ra/((r1 - ra)^2 + (za - zj)^2) + 1)*elliptic_kc(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2)) - elliptic_ec(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2)))*cI0*(2*(r1 - ra)*r1*ra/((r1 - ra)^2 + (za - zj)^2)^2 + r1/((r1 - ra)^2 + (za - zj)^2))/(r1*(4*r1*ra/((r1 - ra)^2 + (za - zj)^2) + 1)*ra) - 2*sqrt((r2 - ra)^2 + (za - zj)^2)*((4*r2*ra/((r2 - ra)^2 + (za - zj)^2) + 1)*elliptic_kc(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2)) - elliptic_ec(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2)))*cI0*(2*(r2 - ra)*r2*ra/((r2 - ra)^2 + (za - zj)^2)^2 + r2/((r2 - ra)^2 + (za - zj)^2))/(r2*(4*r2*ra/((r2 - ra)^2 + (za - zj)^2) + 1)*ra)


file = open('field_of_deyna_cylinder.txt', 'a')
file.write("As_diff_ra(cI0, rj1, rj2, ra, zj, za) = " + str(As_diff_ra(cI0, rj1, rj2, ra, zj, za)))
file.write('\n\n')
file.close()

Iphi_jv_rj_diff_ra = lambda c_I0, r_j, r_a, z_j, z_a : (Iphi_jv_rj.diff(ra)).substitute(ra==r_a).substitute(za==z_a).substitute(zj==z_j).substitute(rj==r_j).substitute(cI0==c_I0)
print "Iphi_jv_rj_diff_ra(cI0, rj, ra, zj, za)  =", Iphi_jv_rj_diff_ra(cI0, rj, ra, zj, za)

file = open('field_of_deyna_cylinder.txt', 'a')
file.write("Iphi_jv_rj_diff_ra(cI0, rj, ra, zj, za) = " + str(Iphi_jv_rj_diff_ra(cI0, rj, ra, zj, za)))
file.write('\n\n')
file.close()

#Av_diff_ra = integrate(Iphi_jv_rj_diff_ra, (rj, rj1, rj2), algorithm="giac")
Av_diff_ra = lambda cI0, rj1, rj2, ra, zj, za : my_numerical_integral(lambda rj : Iphi_jv_rj_diff_ra(cI0, rj, ra, zj, za), rj1, rj2)
# print "Av_diff_ra  =", Av_diff_ra(cI0, rj1, rj2, ra, zj, za)


epsilon = 0.000
Zj1 = -3.0 + epsilon
Zj2 =  3.0 - epsilon
Rj1 = 0.5 + epsilon
Rj2 = 1.0  - epsilon

# sizes of cylinders in deyna video
Zj1 = -1.5
Zj2 =  1.5
Rj1 = 0.3
Rj2 = 1.5

Ra1 = Rj1
Ra2 = Rj2

DZ = Zj2 - Zj1

file = open('field_of_deyna_cylinder.txt', 'a')
file.write("sizes of cylinder")
file.write('\n')
file.write("Zj1 = " + str(Zj1) + ", Zj2 = " + str(Zj2))
file.write('\n')
file.write("Rj1 = " + str(Rj1) + ", Rj1 = " + str(Rj1))
file.write('\n')
file.write("full_volume_cylinder = " + str(full_volume_cylinder))
file.write('\n\n')
file.close()

# At_diff_za_substituted_zj = At_diff_za.substitute(cI0==1, zj1==Zj1, zj2==Zj2)
# At_diff_za_substituted_zj = At_diff_za(cI0, rj, ra, zj1, zj2, za).substitute(cI0==1, zj1==Zj1, zj2==Zj2)
At_diff_za_substituted_zj = lambda rj, ra, za : At_diff_za(1, rj, ra, Zj1, Zj2, za)
# print "At_diff_za_substituted_zj =", At_diff_za_substituted_zj

file = open('field_of_deyna_cylinder.txt', 'a')
file.write('At_diff_za_substituted_zj(rj, ra, za) = ' + str(At_diff_za_substituted_zj(rj, ra, za)))
file.write('\n\n')
file.write('At_diff_za_substituted_zj(1.0, 1.1, 0.9) = ' + str(At_diff_za_substituted_zj(1.0, 1.1, 0.9)))
file.write('\n\n')
file.close()

# As_diff_ra_substituted_rj = As_diff_ra.substitute(cI0==1, rj1==Rj1, rj2==Rj2)
# As_diff_ra_substituted_rj = As_diff_ra(cI0, rj1, rj2, ra, zj, za).substitute(cI0==1, rj1==Rj1, rj2==Rj2)
As_diff_ra_substituted_rj = lambda ra, zj, za : As_diff_ra(1, Rj1, Rj2, ra, zj, za)
# print "As_diff_ra_substituted_rj =", As_diff_ra_substituted_rj

file = open('field_of_deyna_cylinder.txt', 'a')
file.write('As_diff_ra_substituted_rj(ra, zj, za) = ' + str(As_diff_ra_substituted_rj(ra, zj, za)))
file.write('\n\n')
file.write('As_diff_ra_substituted_rj(1.1, 0.1, 0.9) = ' + str(As_diff_ra_substituted_rj(1.1, 0.1, 0.9)))
file.write('\n\n')
file.close()

# Av_diff_ra.substitute(cI0==1, rj1==Rj1, rj2==Rj2)
# Av_diff_ra_substituted_rj = lambda ra, zj, za : Av_diff_ra (cI0, rj1, rj2, ra, zj, za).substitute(cI0==1, rj1==Rj1, rj2==Rj2)
Av_diff_ra_substituted_rj = lambda ra, zj, za : Av_diff_ra (1, Rj1, Rj2, ra, zj, za)
print "Av_diff_ra_substituted_rj(1.1, 0.1, 0.9) =", Av_diff_ra_substituted_rj(1.1, 0.1, 0.9)

file = open('field_of_deyna_cylinder.txt', 'a')
#file.write('Av_diff_ra_substituted_rj(ra, zj, za) = ' + str(Av_diff_ra_substituted_rj(ra, zj, za)))
#file.write('\n\n')
file.write('Av_diff_ra_substituted_rj(1.1, 0.1, 0.9) = ' + str(Av_diff_ra_substituted_rj(1.1, 0.1, 0.9)))
file.write('\n\n')
file.close()

def calc_H_phi( At_diff_za_substituted_zj, As_diff_ra_substituted_rj, Za, Ra):
    # At_diff_za_substituted_zj_za_ra = At_diff_za_substituted_zj.substitute(za==Za, ra==Ra)
    # At_diff_za_substituted_zj_za_ra = At_diff_za_substituted_zj(rj, ra, za).substitute(za==Za, ra==Ra)
    At_diff_za_substituted_zj_za_ra = lambda rj : At_diff_za_substituted_zj(rj, Ra, Za)
    # print "At_diff_za_substituted_zj_za_ra =", At_diff_za_substituted_zj_za_ra

    # As_diff_ra_substituted_rj_za_ra = As_diff_ra_substituted_rj.substitute(za==Za, ra==Ra)
    # As_diff_ra_substituted_rj_za_ra = As_diff_ra_substituted_rj(ra, zj, za).substitute(za==Za, ra==Ra)
    As_diff_ra_substituted_rj_za_ra = lambda zj : As_diff_ra_substituted_rj(Ra, zj, Za)
    # print "As_diff_ra_substituted_rj_za_ra =", As_diff_ra_substituted_rj_za_ra

    # Av_diff_ra_substituted_rj_za_ra = Av_diff_ra_substituted_rj.substitute(za==Za, ra==Ra)
    # Av_diff_ra_substituted_rj_za_ra = lambda zj : Av_diff_ra_substituted_rj(Ra, zj, Za)
    # print "Av_diff_ra_substituted_rj_za_ra =", Av_diff_ra_substituted_rj_za_ra

    # debug
    print At_diff_za_substituted_zj_za_ra(rj).substitute(rj == Rj2)
    print As_diff_ra_substituted_rj_za_ra(zj).substitute(zj == 0)
    # print Av_diff_ra_substituted_rj_za_ra.substitute(zj == 0)
    # print Av_diff_ra_substituted_rj_za_ra(0)
    print Av_diff_ra_substituted_rj(Ra, 0, Za)

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write("debug calc_H_phi")
    file.write('\n')
    file.write("At_diff_za_substituted_zj_za_ra(rj).substitute(rj == Rj2) = " + str(At_diff_za_substituted_zj_za_ra(rj).substitute(rj == Rj2)))
    file.write('\n')
    file.write("As_diff_ra_substituted_rj_za_ra(zj).substitute(zj == 0) = " + str(As_diff_ra_substituted_rj_za_ra(zj).substitute(zj == 0)))
    file.write('\n')
    file.write("Av_diff_ra_substituted_rj(Ra, 0, Za) = " + str(Av_diff_ra_substituted_rj(Ra, 0, Za)))
    file.write('\n\n')
    file.close()

    At_diff_za_num_int = At_diff_za_substituted_zj_za_ra(rj).nintegral(rj, Rj1, Rj2)
    # print "At_diff_za_num_int  =", At_diff_za_num_int

    # As_v_diff_ra_num_int = (As_diff_ra_substituted_rj_za_ra - Av_diff_ra_substituted_rj_za_ra  ).nintegral(zj, Zj1, Zj2)
    As_diff_ra_num_int = (As_diff_ra_substituted_rj_za_ra(zj)).nintegral(zj, Zj1, Zj2)

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write("Av_diff_ra_num_int = my_numerical_integral( lambda zj : Av_diff_ra_substituted_rj(Ra, zj, Za), Zj1, Zj2)")
    file.write('\n\n')
    file.close()

    # Av_diff_ra_num_int = my_numerical_integral( lambda zj : Av_diff_ra_substituted_rj_za_ra(zj), Zj1, Zj2)
    Av_diff_ra_num_int = my_numerical_integral( lambda zj : Av_diff_ra_substituted_rj(Ra, zj, Za), Zj1, Zj2)
    print "Av_diff_ra_num_int =", Av_diff_ra_num_int

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write("Av_diff_ra_num_int = " + str(Av_diff_ra_num_int))
    file.write('\n\n')
    file.close()

    As_v_diff_ra_num_int = As_diff_ra_num_int[0] - Av_diff_ra_num_int
    # print "As_v_diff_ra_num_int  =", As_v_diff_ra_num_int

    H_phi_t = At_diff_za_num_int[0]
    H_phi_s = - As_v_diff_ra_num_int

    H_phi = H_phi_t + H_phi_s
    print "Ra=", Ra, "Za=", Za, "H_phi_t =", H_phi_t
    print "Ra=", Ra, "Za=", Za, "H_phi_s =", H_phi_s
    print "Ra=", Ra, "Za=", Za, "H_phi =", H_phi

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write('\n')
    file.write("Ra=" + str(Ra) + ", Za=" + str(Za) + ", H_phi_t=" + str(H_phi_t))
    file.write('\n')
    file.write("Ra=" + str(Ra) + ", Za=" + str(Za) + ", H_phi_s=" + str(H_phi_s))
    file.write('\n')
    file.write("Ra=" + str(Ra) + ", Za=" + str(Za) + ", H_phi  =" + str(H_phi))
    file.write('\n\n')
    file.close()

    return (H_phi, H_phi_t, H_phi_s)


def calc_F_lorenz( At_diff_za_substituted_zj, As_diff_ra_substituted_rj, Za, Ra1, Ra2):
    # At_diff_za_substituted_zj_za = lambda rj, ra : At_diff_za_substituted_zj.substitute(za==Za)
    # At_diff_za_substituted_zj_za = lambda rj, ra : At_diff_za_substituted_zj(rj, ra, za).substitute(za==Za)
    At_diff_za_substituted_zj_za = lambda rj, ra : At_diff_za_substituted_zj(rj, ra, Za)
    print "At_diff_za_substituted_zj_za =", At_diff_za_substituted_zj_za

    # As_diff_ra_substituted_rj_za = lambda ra, zj : As_diff_ra_substituted_rj.substitute(za==Za)
    # As_diff_ra_substituted_rj_za = lambda ra, zj : As_diff_ra_substituted_rj(ra, zj, za).substitute(za==Za)
    As_diff_ra_substituted_rj_za = lambda ra, zj : As_diff_ra_substituted_rj(ra, zj, Za)
    print "As_diff_ra_substituted_rj_za =", As_diff_ra_substituted_rj_za

    # Av_diff_ra_substituted_rj_za = Av_diff_ra_substituted_rj.substitute(za==Za)
    Av_diff_ra_substituted_rj_za = lambda ra, zj : Av_diff_ra_substituted_rj(ra, zj, Za)
    # print "Av_diff_ra_substituted_rj_za =", Av_diff_ra_substituted_rj_za

    # debug
    print At_diff_za_substituted_zj_za(rj, ra)
    print As_diff_ra_substituted_rj_za(ra, zj)

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write("debug calc_F_lorenz")
    file.write('\n\n')
    file.write("At_diff_za_substituted_zj_za(rj, ra) = " + str(At_diff_za_substituted_zj_za(rj, ra)))
    file.write('\n\n')
    file.write("As_diff_ra_substituted_rj_za(ra, zj) = " + str(As_diff_ra_substituted_rj_za(ra, zj)))
    file.write('\n\n')
    file.write("Av_diff_ra_substituted_rj_za((Ra1 + Ra2)/2, 0) = " + str(Av_diff_ra_substituted_rj_za((Ra1 + Ra2)/2, 0)))
    file.write('\n\n')
    file.close()

    jt_substituted_cI0 = jt.substitute(cI0 == 1)

    At_diff_za_num_int_ra = lambda Rj : my_numerical_integral(lambda ra : (2*pi*jt_substituted_cI0*ra*At_diff_za_substituted_zj_za(Rj, ra)), Ra1, Ra2)

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write("At_diff_za_num_int_ra_int_rj = my_numerical_integral(lambda rj : At_diff_za_num_int_ra(rj), Rj1, Rj2)")
    # file.write("At_diff_za_num_int_ra(rj) = " + str(At_diff_za_num_int_ra(rj)))
    file.write('\n\n')
    file.close()

    At_diff_za_num_int_ra_int_rj = my_numerical_integral(lambda rj : At_diff_za_num_int_ra(rj), Rj1, Rj2)
    print "At_diff_za_num_int_ra_int_rj =", At_diff_za_num_int_ra_int_rj

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write("At_diff_za_num_int_ra_int_rj = " + str(At_diff_za_num_int_ra_int_rj))
    file.write('\n\n')
    file.close()

    #As_v_diff_ra_num_int_int = (As_diff_ra_substituted_rj_za - Av_diff_ra  ).nintegral(ra,Ra1,Ra2).nintegral(zj, Zj1, Zj2)
    As_diff_ra_num_int_ra = lambda Zj : my_numerical_integral(lambda ra : (2*pi*jt_substituted_cI0*ra*(As_diff_ra_substituted_rj_za(ra, Zj) ) ), Ra1, Ra2)

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write("As_diff_ra_num_int_ra_int_zj = my_numerical_integral(lambda zj : As_diff_ra_num_int_ra(zj), Zj1, Zj2)")
    #file.write("As_diff_ra_num_int_ra(zj) = " + str(As_diff_ra_num_int_ra(zj)))
    file.write('\n\n')
    file.close()

    As_diff_ra_num_int_ra_int_zj = my_numerical_integral(lambda zj : As_diff_ra_num_int_ra(zj), Zj1, Zj2)
    print "As_diff_ra_num_int_ra_int_zj  =", As_diff_ra_num_int_ra_int_zj

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write("As_diff_ra_num_int_ra_int_zj = " + str(As_diff_ra_num_int_ra_int_zj))
    file.write('\n\n')
    file.close()

    Av_diff_ra_num_int_ra = lambda Zj : my_numerical_integral(lambda ra : (2*pi*jt_substituted_cI0*ra*(Av_diff_ra_substituted_rj_za(ra, Zj) ) ), Ra1, Ra2)

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write("Av_diff_ra_num_int_ra_int_zj = my_numerical_integral(lambda zj : Av_diff_ra_num_int_ra(zj), Zj1, Zj2)")
    file.write('\n\n')
    file.close()

    Av_diff_ra_num_int_ra_int_zj = my_numerical_integral(lambda zj : Av_diff_ra_num_int_ra(zj), Zj1, Zj2)
    print "Av_diff_ra_num_int_ra_int_zj  =", Av_diff_ra_num_int_ra_int_zj

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write("Av_diff_ra_num_int_ra_int_zj = " + str(Av_diff_ra_num_int_ra_int_zj))
    file.write('\n\n')
    file.close()

    F_z_t = At_diff_za_num_int_ra_int_rj
    F_z_s = As_diff_ra_num_int_ra_int_zj
    F_z_v = Av_diff_ra_num_int_ra_int_zj

    As_v_diff_ra_num_int_ra_int_zj = As_diff_ra_num_int_ra_int_zj - Av_diff_ra_num_int_ra_int_zj
    F_z_sv = - As_v_diff_ra_num_int_ra_int_zj
    F_z = F_z_t + F_z_sv

    print "Ra1=", Ra1, "Ra2=", Ra2, "Za=", Za, "F_z_t =", F_z_t
    print "Ra1=", Ra1, "Ra2=", Ra2, "Za=", Za, "F_z_s =", F_z_s
    print "Ra1=", Ra1, "Ra2=", Ra2, "Za=", Za, "F_z_v =", F_z_v
    print "Ra1=", Ra1, "Ra2=", Ra2, "Za=", Za, "F_z_sv =", F_z_sv
    print "Ra1=", Ra1, "Ra2=", Ra2, "Za=", Za, "F_z =", F_z

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write('\n')
    file.write("Ra1=" + str(Ra1) + ", Ra2=" + str(Ra2) + ", Za=" + str(Za) + ", F_z_t = " + str(F_z_t))
    file.write('\n')
    file.write("Ra1=" + str(Ra1) + ", Ra2=" + str(Ra2) + ", Za=" + str(Za) + ", F_z_s = " + str(F_z_s))
    file.write('\n')
    file.write("Ra1=" + str(Ra1) + ", Ra2=" + str(Ra2) + ", Za=" + str(Za) + ", F_z_v = " + str(F_z_v))
    file.write('\n')
    file.write("Ra1=" + str(Ra1) + ", Ra2=" + str(Ra2) + ", Za=" + str(Za) + ", F_z_sv = " + str(F_z_sv))
    file.write('\n')
    file.write("Ra1=" + str(Ra1) + ", Ra2=" + str(Ra2) + ", Za=" + str(Za) + ", F_z   = " + str(F_z))
    file.write('\n\n')
    file.close()

    return (F_z, F_z_t, F_z_s, F_z_v)

def calc_F_lorenz_cylinder(dz):
    # расчет силы Лоренца, действующей на ближайжий (правый) торец пробного цилиндра расположенного левее на расстоянии
    Za = Zj1 - dz
    F_lorenz_left_cylinder_right_t = calc_F_lorenz( At_diff_za_substituted_zj, As_diff_ra_substituted_rj, Za, Ra1, Ra2)
    print "F_lorenz_left_cylinder_right_t =", F_lorenz_left_cylinder_right_t

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write("calc_F_lorenz_cylinder dz = " + str(dz) + ", Za = " + str(Za))
    file.write('\n')
    file.write("F_lorenz_left_cylinder_right_t = " + str(F_lorenz_left_cylinder_right_t))
    file.write('\n\n')
    file.close()

    # расчет силы Лоренца, действующей на удалённый (левый) торец пробного цилиндра расположенного левее на расстоянии
    Za = Zj1 - DZ - dz
    F_lorenz_left_cylinder_left_t = calc_F_lorenz( At_diff_za_substituted_zj, As_diff_ra_substituted_rj, Za, Ra1, Ra2)
    print "F_lorenz_left_cylinder_left_t =", F_lorenz_left_cylinder_left_t

    F_lorenz_cylinder = - F_lorenz_left_cylinder_right_t[0] + F_lorenz_left_cylinder_left_t[0]
    print "F_lorenz_cylinder =", F_lorenz_cylinder

    file = open('field_of_deyna_cylinder.txt', 'a')
    file.write("calc_F_lorenz_cylinder dz = " + str(dz) + ", Za = " + str(Za))
    file.write('\n')
    file.write("F_lorenz_left_cylinder_left_t = " + str(F_lorenz_left_cylinder_left_t))
    file.write('\n')
    file.write("F_lorenz_cylinder = " + str(F_lorenz_cylinder))
    file.write('\n\n')
    file.close()

    return F_lorenz_cylinder

dr = 0.01
dz = 0.04

'''
Za = 0
Ra = (Rj1 + Rj2) / 2
calc_H_phi( At_diff_za_substituted_zj, As_diff_ra_substituted_rj, Za, Ra)

Za = 0
Ra = Rj2 - dr
calc_H_phi( At_diff_za_substituted_zj, As_diff_ra_substituted_rj, Za, Ra)

Za = 0
Ra = Rj1 + dr
calc_H_phi( At_diff_za_substituted_zj, As_diff_ra_substituted_rj, Za, Ra)
'''

plot_data_f = []
plot_data_h = []
plot_data_h_t = []
plot_data_h_s = []

Ra = (Rj1 + Rj2) / 2
for dz in (0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0):
    Za = Zj1 - dz
    h = calc_H_phi( At_diff_za_substituted_zj, As_diff_ra_substituted_rj, Za, Ra)
    f = calc_F_lorenz_cylinder(dz)
    plot_data_h += [(Za, h[0])]
    plot_data_h_t += [(Za, h[1])]
    plot_data_h_s += [(Za, h[2])]
    plot_data_f += [(Za, f)]


dir = os.getcwd()  + "/results/horizontal_test_plot"
print "dir = ", dir

try:
    os.mkdir(dir)
except:
    pass

folder = dir + "/"

params = \
    "_Rj1=" + float_formatting(Rj1) + \
    "_Rj2=" + float_formatting(Rj2) + \
    "_Zj1=" + float_formatting(Zj1) + \
    "_Zj2=" + float_formatting(Zj2) + \
    "_full_volume=" + str(full_volume_cylinder) + \
    "_Ra=" + float_formatting(Ra)

p = list_plot(plot_data_f)
pname = folder + "F_lorenz_cylinder" + params + ".png"
print pname
p.save(pname)

p = list_plot(plot_data_h)
pname = folder + "H_phi" + params + ".png"
print pname
p.save(pname)

p = list_plot(plot_data_h_t)
pname = folder + "H_phi_t" + params + ".png"
print pname
p.save(pname)

p = list_plot(plot_data_h_s)
pname = folder + "H_phi_s" + params + ".png"
print pname
p.save(pname)

"""

epsilon = 0.000
Zj1 = -3.0 + epsilon
Zj2 =  3.0 - epsilon
Rj1 = 0.5 + epsilon
Rj2 = 1.0  - epsilon

-1.11938314635325
23.5205877480921
Ra= 0.750000000000000 Za= 0 H_phi_t = -0.590667166965
Ra= 0.750000000000000 Za= 0 H_phi_s = -16.7063519209At_diff_za_substituted
Ra= 0.750000000000000 Za= 0 H_phi = -17.2970190879
-1.06877972270768
203.282584697895
Ra= 0.990000000000000 Za= 0 H_phi_t = -0.559994640818
Ra= 0.990000000000000 Za= 0 H_phi_s = -12.6345175482
Ra= 0.990000000000000 Za= 0 H_phi = -13.194512189
-1.15783939892903
417.808854381580
Ra= 0.510000000000000 Za= 0 H_phi_t = -0.614427172046
Ra= 0.510000000000000 Za= 0 H_phi_s = -24.6044152761
Ra= 0.510000000000000 Za= 0 H_phi = -25.2188424482
0.212265859504912
0.0278735129943648
Ra= 0.750000000000000 Za= -3.01000000000000 H_phi_t = 8.09788018021
Ra= 0.750000000000000 Za= -3.01000000000000 H_phi_s = -8.14062009287
Ra= 0.750000000000000 Za= -3.01000000000000 H_phi = -0.0427399126641
0.584992775274781
0.0274843166541227
Ra= 0.750000000000000 Za= -3.02000000000000 H_phi_t = 7.90079601518
Ra= 0.750000000000000 Za= -3.02000000000000 H_phi_s = -7.90627591381
Ra= 0.750000000000000 Za= -3.02000000000000 H_phi = -0.00547989862735
0.950751444665585
0.0271014637703293
Ra= 0.750000000000000 Za= -3.03000000000000 H_phi_t = 7.70374457824
Ra= 0.750000000000000 Za= -3.03000000000000 H_phi_s = -7.67337390159
Ra= 0.750000000000000 Za= -3.03000000000000 H_phi = 0.0303706766451
1.30632991706071
0.0267248369795651
Ra= 0.750000000000000 Za= -3.04000000000000 H_phi_t = 7.50745634914
Ra= 0.750000000000000 Za= -3.04000000000000 H_phi_s = -7.44260170471
Ra= 0.750000000000000 Za= -3.04000000000000 H_phi = 0.0648546444295
1.64882183140010
0.0263543213333033
Ra= 0.750000000000000 Za= -3.05000000000000 H_phi_t = 7.31262434998
Ra= 0.750000000000000 Za= -3.05000000000000 H_phi_s = -7.21460852535
Ra= 0.750000000000000 Za= -3.05000000000000 H_phi = 0.0980158246291
3.09180381421842
0.0245895493685892_data_h = []
plot_data_h_t = []
plot_d
Ra= 0.750000000000000 Za= -3.10000000000000 H_phi_t = 6.38079682334
Ra= 0.750000000000000 Za= -3.10000000000000 H_phi_s = -6.13526547968
Ra= 0.750000000000000 Za= -3.10000000000000 H_phi = 0.24553134366
4.46844702245850
0.0214585502234302
Ra= 0.750000000000000 Za= -3.20000000000000 H_phi_t = 4.85596237851
Ra= 0.750000000000000 Za= -3.20000000000000 H_phi_s = -4.3904039877
Ra= 0.750000000000000 Za= -3.20000000000000 H_phi = 0.465558390813
4.54322990976506
0.0187848913224443
Ra= 0.750000000000000 Za= -3.30000000000000 H_phi_t = 3.7975290884
Ra= 0.750000000000000 Za= -3.300000000000   H_phi_s = -3.18494481861
Ra= 0.750000000000000 Za= -3.30000000000000 H_phi = 0.612584269794
4.17416247979194
0.0164939225375962
Ra= 0.750000000000000 Za= -3.40000000000000 H_phi_t = 3.0727722543
Ra= 0.750000000000000 Za= -3.40000000000000 H_phi_s = -2.3646207498
Ra= 0.750000000000000 Za= -3.40000000000000 H_phi = 0.708151504504
3.72746212762838
0.0145243546693854
Ra= 0.750000000000000 Za= -3.50000000000000 H_phi_t = 2.55928357008
Ra= 0.750000000000000 Za= -3.50000000000000 H_phi_s = -1.79335064128
Ra= 0.750000000000000 Za= -3.50000000000000 H_phi = 0.765932928794
2.15308921458995
0.00800288860211396
Ra= 0.750000000000000 Za= -4.00000000000000 H_phi_t = 1.30496473098
Ra= 0.750000000000000 Za= -4.00000000000000 H_phi_s = -0.552627497038
Ra= 0.750000000000000 Za= -4.00000000000000 H_phi = 0.752337233944
1.37561504304263
0.00467298952333006
Ra= 0.750000000000000 Za= -4.50000000000000 H_phi_t = 0.790586086289
Ra= 0.750000000000000 Za= -4.50000000000000 H_phi_s = -0.20939598622
Ra= 0.750000000000000 Za= -4.50000000000000 H_phi = 0.581190100069
0.933472562092955
0.00286443942825520
Ra= 0.750000000000000 Za= -5.00000000000000 H_phi_t = 0.517524861151
Ra= 0.750000000000000 Za= -5.00000000000000 H_phi_s = -0.0915423141092
Ra= 0.750000000000000 Za= -5.00000000000000 H_phi = 0.425982547041
1.26118266458605
0.0114998224453978
Ra= 1.50000000000000 Za= -4.00000000000000 H_phi_t = 0.574067356272
Ra= 1.50000000000000 Za= -4.00000000000000 H_phi_s = -0.295288120362
Ra= 1.50000000000000 Za= -4.00000000000000 H_phi = 0.278779235911
"""

#At_diff_za_num_int  = (0.892985470996911, 9.91413030556511e-15, 21, 0)
#As_v_diff_ra_num_int  = (0.7922729550346237, 1.915043453188808e-14, 21, 0)
#H_phi  = 0.100712515962



step_Za = 0.01
min_Za = Zj1-2
max_Za = Zj1


#for Ra in (Rj1/2, Rj1-dr, Rj1, (Rj1 + Rj2) / 2, Rj2, Rj2+dr, Rj2+Rj1/2, Rj2+(Rj1 + Rj2) / 2):
#for Ra in ( Rj2, Rj2+dr, Rj2+Rj1/2, Rj2+(Rj1 + Rj2) / 2):
for Ra in np.arange(Rj1, Rj2+0.1, 0.1):
    plot_data_h = []
    plot_data_h_t = []
    plot_data_h_s = []
    for Za in np.arange(min_Za, max_Za, step_Za):
        try:
            h = calc_H_phi( At_diff_za_substituted_zj, As_diff_ra_substituted_rj, Za, Ra)
            plot_data_h += [(Za, h[0])]
            plot_data_h_t += [(Za, h[1])]
            plot_data_h_s += [(Za, h[2])]
        except:
            pass

    dir = os.getcwd()  + "/results/horizontal_plot"
    print "dir = ", dir

    try:
        os.mkdir(dir)
    except:
        pass

    folder = dir + "/"

    params = \
        "_Rj1=" + float_formatting(Rj1) + \
        "_Rj2=" + float_formatting(Rj2) + \
        "_Zj1=" + float_formatting(Zj1) + \
        "_Zj2=" + float_formatting(Zj2) + \
        "_full_volume=" + str(full_volume_cylinder) + \
        "_Ra=" + float_formatting(Ra) + \
        "_min_Za=" + float_formatting(min_Za) + \
        "_max_Za=" + float_formatting(max_Za)

    p = list_plot(plot_data_h)
    pname = folder + "H_phi" + params + ".png"
    print pname
    p.save(pname)

    p = list_plot(plot_data_h_t)
    pname = folder + "H_phi_t" + params + ".png"
    print pname
    p.save(pname)

    p = list_plot(plot_data_h_s)
    pname = folder + "H_phi_s" + params + ".png"
    print pname
    p.save(pname)


'''



step_Ra = 0.1
min_Ra = 0.1
max_Ra = 1.5

ra_range = [min_Ra, max_Ra]
ra_range2 = [-max_Ra, max_Ra]
n_ra = round((max_Ra-min_Ra)/step_Ra) + 1
print "n_ra =", n_ra

step_Za = 0.1
min_Za = -5
max_Za = 5

za_range = [min_Za, max_Za]
n_za = round((max_Za-min_Za)/step_Za) + 1
print "n_za =", n_za


plot_data = []
arr_H = np.zeros((n_ra*2+1, n_za))
nr = 0
for Ra in np.arange(min_Ra, max_Ra + step_Ra, step_Ra):
    n_nan = 0
    h = []
    nz = 0
    for Za in np.arange(min_Za, max_Za + step_Za, step_Za):
        try:
            At_diff_za_substituted2 = At_diff_za_substituted_zj.substitute(za==Za, ra==Ra)
            #print "At_diff_za_substituted2 =", At_diff_za_substituted2
            As_diff_ra_substituted2 = As_diff_ra_substituted_rj.substitute(za==Za, ra==Ra)
            #print "As_diff_ra_substituted2 =", As_diff_ra_substituted2

            At_diff_za_num_int = At_diff_za_substituted2.nintegral(rj, R1, R2)
            #print "At_diff_za_num_int  =", At_diff_za_num_int

            As_v_diff_ra_num_int = (As_diff_ra_substituted2 - Av_diff_ra  ).nintegral(zj, Z1, Z2)
            #print "As_v_diff_ra_num_int  =", As_v_diff_ra_num_int

            H_phi = At_diff_za_num_int[0] - As_v_diff_ra_num_int[0]
            print "Ra  =", Ra, "Za  =", Za, "H_phi  =", H_phi

            plot_data.append((Za, Ra, H_phi))
            plot_data.append((Za, -Ra, -H_phi))
            h.append(H_phi)
        except:
            H_phi = np.nan
            n_nan = n_nan + 1
        print "nr=", nr
        print "nz=", nz
        arr_H[nr,nz] = H_phi
        arr_H[n_ra*2 - nr,nz] = - H_phi
        nz = nz + 1
    nr = nr + 1
    print "len(h)=", len(h)
    print "n_nan=", n_nan

print "nr=", nr
print "nz=", nz

results_folder = "../articles/sagemath/field_of_deyna_cylinder/results/"
results_folder = "./results/"

print plot_data
print arr_H

from sage.plot.contour_plot import ContourPlot
p = ContourPlot(arr_H, za_range, ra_range2, options=dict(fill=True, contours=np.arange(-10, 10, 0.5)))
g = Graphics()
g.add_primitive(p)

from sage.plot.line import Line
from sage.plot.arrow import Arrow

L = Line([-5,5],[0,0],{'alpha':1,'thickness':1,'rgbcolor':(0,1,1),'legend_label':''})
g.add_primitive(L)

A = Arrow(Z2,R1,Z1,R1,{'width':1,'head':1,'rgbcolor':(1,0,0),'linestyle':'dashed','zorder':8,'legend_label':''})
g.add_primitive(A)

A = Arrow(Z1,R2,Z2,R2,{'width':1,'head':1,'rgbcolor':(1,0,0),'linestyle':'dashed','zorder':8,'legend_label':''})
g.add_primitive(A)

A = Arrow(Z1,R1,Z1,R2,{'width':1,'head':1,'rgbcolor':(1,0,0),'linestyle':'dashed','zorder':8,'legend_label':''})
g.add_primitive(A)

A = Arrow(Z2,R2,Z2,R1,{'width':1,'head':1,'rgbcolor':(1,0,0),'linestyle':'dashed','zorder':8,'legend_label':''})
g.add_primitive(A)



pname = results_folder + "H_phi_contour.09" + ".png"
print pname
g.save(pname)


g = list_plot3d(plot_data, frame_aspect_ratio=[1, 1, 1/3])


pname = results_folder + "H_phi.09" + ".png"
print pname
g.save(pname)

for d in plot_data:
    print d

for d in plot_data:
    if d[1] > 0.5 and d[1] < 0.7:
        print "za=",d[0],"ra=", d[1], "H_phi=", d[2]
# '''
















