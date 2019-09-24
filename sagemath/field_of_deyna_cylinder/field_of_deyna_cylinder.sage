import sys
reload(sys)
sys.setdefaultencoding('utf8')


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



print elliptic_kc(0.5)
# 1.85407467730137

print elliptic_kc(-0.5)
# 1.41573720842596

m = var("m")
print elliptic_kc(m).diff(m)
# -1/2*((m - 1)*elliptic_kc(m) + elliptic_ec(m))/((m - 1)*m)

rja2 = (rj-ra)^2+(zj-za)^2
module = - 4*rj*ra / rja2
Iphi=4*elliptic_kc(module)/sqrt(rja2)
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


# строго говоря размерность поверхностного тока должна быть приведена в соответствие с формулой js = c*[I x n]
# а размерность объёмного - формуле jv = c * rot(I)

cI0 = var("cI0")
js = cI0 / rj
jv = 0

#js = cI0
#jv = cI0 / rj

z1 = var("z1")
z2 = var("z2")

r1 = var("r1")
r2 = var("r2")


Iphi_js_rj=Iphi*js*rj
print "Iphi_js_rj =", Iphi_js_rj

Iphi_jv_rj=Iphi*jv*rj
print "Iphi_jv_rj =", Iphi_jv_rj

# As = Iphi_js_rj.integrate(zj, algorithm="mathematica_free")
# NotImplementedError: Mathematica online integrator can only handle single letter variables.


# As = integrate(Iphi_js_rj, zj, algorithm="sympy")
# infinity loop

# As = integrate(Iphi_js_rj, zj, algorithm="fricas")
# TypeError: unable to start fricas because the command 'fricas -nosman' failed: The command was not found or was not executable: fricas.

# maple(Iphi_js_rj).integrate(zj)

# As = integrate(Iphi_js_rj, zj)
# RuntimeError: Encountered operator mismatch in maxima-to-sr translation

# As = integrate(Iphi_js_rj, zj, algorithm="giac")

As = integrate(Iphi_js_rj, (zj, z1, z2), algorithm="giac")
print "As =", As
# As = integrate(4*rj*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (za - zj)^2))/sqrt((ra - rj)^2 + (za - zj)^2), zj)

At = integrate(Iphi_js_rj, (rj, r1, r2), algorithm="giac")
print "At =", At
# At = integrate(4*cI0*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (za - zj)^2))/sqrt((ra - rj)^2 + (za - zj)^2), rj, r1, r2)

Av = integrate(integrate(Iphi_jv_rj, (zj, z1, z2), algorithm="giac"), (rj, r1, r2), algorithm="giac")
print "Av =", Av

# AT = At.substitute_expression(zj==z1) - At.substitute_expression(zj==z2)
# RuntimeError: Encountered operator mismatch in maxima-to-sr translation

AT = integrate(Iphi_js_rj.substitute_expression(zj==z1), (rj, r1, r2), algorithm="giac") - integrate(Iphi_js_rj.substitute_expression(zj==z2), (rj, r1, r2), algorithm="giac")
print "AT =", AT










