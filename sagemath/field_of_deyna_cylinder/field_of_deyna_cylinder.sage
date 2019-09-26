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

# dIphi_dza = Iphi.diff(za)
# print "dIphi_dza =", dIphi_dza
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

# dIphi_dra = Iphi.diff(ra)
# print "dIphi_dra =", dIphi_dra
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

# As = integrate(Iphi_js_rj, (zj, z1, z2), algorithm="giac")
# print "As =", As
# As = integrate(4*rj*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (za - zj)^2))/sqrt((ra - rj)^2 + (za - zj)^2), zj)

# At = integrate(Iphi_js_rj, (rj, r1, r2), algorithm="giac")
# print "At =", At
# At = integrate(4*cI0*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (za - zj)^2))/sqrt((ra - rj)^2 + (za - zj)^2), rj, r1, r2)

# Av = integrate(integrate(Iphi_jv_rj, (zj, z1, z2), algorithm="giac"), (rj, r1, r2), algorithm="giac")
# print "Av =", Av
# Av = 0

# AT = At.substitute_expression(zj==z1) - At.substitute_expression(zj==z2)
# RuntimeError: Encountered operator mismatch in maxima-to-sr translation

# AT = integrate(Iphi_js_rj.substitute(zj==z1), (rj, r1, r2), algorithm="giac") - integrate(Iphi_js_rj.substitute(zj==z2), (rj, r1, r2), algorithm="giac")
# print "AT =", AT
# AT = integrate(4*cI0*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2))/sqrt((ra - rj)^2 + (z1 - za)^2), rj, r1, r2) - integrate(4*cI0*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2))/sqrt((ra - rj)^2 + (z2 - za)^2), rj, r1, r2)

# AS = - As.substitute_expression(rj==r1) + As.substitute_expression(rj==r2)
# AS =  - integrate(Iphi_js_rj.substitute(rj==r1), (zj, z1, z2), algorithm="giac") + integrate(Iphi_js_rj.substitute(rj==r2), (zj, z1, z2), algorithm="giac")
# print "AS =", AS
# AS = -integrate(4*cI0*elliptic_kc(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2))/sqrt((r1 - ra)^2 + (za - zj)^2), zj, z1, z2) + integrate(4*cI0*elliptic_kc(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2))/sqrt((r2 - ra)^2 + (za - zj)^2), zj, z1, z2)

# H_phi = AT.diff(za) - (AS - Av).diff(ra)
# RuntimeError: Encountered operator mismatch in maxima-to-sr translation

At_diff_za = Iphi_js_rj.diff(za) .substitute(zj==z1) - Iphi_js_rj.diff(za) .substitute(zj==z2)
print "At_diff_za =", At_diff_za
# At_diff_za = 4*cI0*(z1 - za)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2))/((ra - rj)^2 + (z1 - za)^2)^(3/2) - 4*cI0*(z2 - za)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2))/((ra - rj)^2 + (z2 - za)^2)^(3/2) - 4*((4*ra*rj/((ra - rj)^2 + (z1 - za)^2) + 1)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2)) - elliptic_ec(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2)))*cI0*(z1 - za)/(((ra - rj)^2 + (z1 - za)^2)^(3/2)*(4*ra*rj/((ra - rj)^2 + (z1 - za)^2) + 1)) + 4*((4*ra*rj/((ra - rj)^2 + (z2 - za)^2) + 1)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2)) - elliptic_ec(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2)))*cI0*(z2 - za)/(((ra - rj)^2 + (z2 - za)^2)^(3/2)*(4*ra*rj/((ra - rj)^2 + (z2 - za)^2) + 1))

# AT_diff_za = integrate(At_diff_za, (rj, r1, r2), algorithm="giac")
# print "AT_diff_za =", AT_diff_za
# AT_diff_za = integrate(4*cI0*(z1 - za)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2))/((ra - rj)^2 + (z1 - za)^2)^(3/2) - 4*cI0*(z2 - za)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2))/((ra - rj)^2 + (z2 - za)^2)^(3/2) - 4*((4*ra*rj/((ra - rj)^2 + (z1 - za)^2) + 1)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2)) - elliptic_ec(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2)))*cI0*(z1 - za)/(((ra - rj)^2 + (z1 - za)^2)^(3/2)*(4*ra*rj/((ra - rj)^2 + (z1 - za)^2) + 1)) + 4*((4*ra*rj/((ra - rj)^2 + (z2 - za)^2) + 1)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2)) - elliptic_ec(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2)))*cI0*(z2 - za)/(((ra - rj)^2 + (z2 - za)^2)^(3/2)*(4*ra*rj/((ra - rj)^2 + (z2 - za)^2) + 1)), rj, r1, r2)

# AT_diff_za = integrate(Iphi_js_rj.diff(za) .substitute(zj==z1), (rj, r1, r2), algorithm="giac") - integrate(Iphi_js_rj.diff(za) .substitute(zj==z2), (rj, r1, r2), algorithm="giac")
# AT_diff_za = integrate(4*cI0*(z1 - za)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2))/((ra - rj)^2 + (z1 - za)^2)^(3/2) - 4*((4*ra*rj/((ra - rj)^2 + (z1 - za)^2) + 1)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2)) - elliptic_ec(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2)))*cI0*(z1 - za)/(((ra - rj)^2 + (z1 - za)^2)^(3/2)*(4*ra*rj/((ra - rj)^2 + (z1 - za)^2) + 1)), rj, r1, r2) - integrate(4*cI0*(z2 - za)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2))/((ra - rj)^2 + (z2 - za)^2)^(3/2) - 4*((4*ra*rj/((ra - rj)^2 + (z2 - za)^2) + 1)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2)) - elliptic_ec(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2)))*cI0*(z2 - za)/(((ra - rj)^2 + (z2 - za)^2)^(3/2)*(4*ra*rj/((ra - rj)^2 + (z2 - za)^2) + 1)), rj, r1, r2)

As_diff_ra  =  - Iphi_js_rj.diff(ra).substitute(rj==r1) + Iphi_js_rj.diff(ra).substitute(rj==r2)
print "As_diff_ra  =", As_diff_ra
# As_diff_ra  = -4*cI0*(r1 - ra)*elliptic_kc(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2))/((r1 - ra)^2 + (za - zj)^2)^(3/2) + 4*cI0*(r2 - ra)*elliptic_kc(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2))/((r2 - ra)^2 + (za - zj)^2)^(3/2) + 2*sqrt((r1 - ra)^2 + (za - zj)^2)*((4*r1*ra/((r1 - ra)^2 + (za - zj)^2) + 1)*elliptic_kc(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2)) - elliptic_ec(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2)))*cI0*(2*(r1 - ra)*r1*ra/((r1 - ra)^2 + (za - zj)^2)^2 + r1/((r1 - ra)^2 + (za - zj)^2))/(r1*(4*r1*ra/((r1 - ra)^2 + (za - zj)^2) + 1)*ra) - 2*sqrt((r2 - ra)^2 + (za - zj)^2)*((4*r2*ra/((r2 - ra)^2 + (za - zj)^2) + 1)*elliptic_kc(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2)) - elliptic_ec(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2)))*cI0*(2*(r2 - ra)*r2*ra/((r2 - ra)^2 + (za - zj)^2)^2 + r2/((r2 - ra)^2 + (za - zj)^2))/(r2*(4*r2*ra/((r2 - ra)^2 + (za - zj)^2) + 1)*ra)

# AS_diff_ra  =  integrate(As_diff_ra, (zj, z1, z2), algorithm="giac")
# print "AS_diff_ra  =", AS_diff_ra
# AS_diff_ra  = integrate(-4*cI0*(r1 - ra)*elliptic_kc(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2))/((r1 - ra)^2 + (za - zj)^2)^(3/2) + 4*cI0*(r2 - ra)*elliptic_kc(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2))/((r2 - ra)^2 + (za - zj)^2)^(3/2) + 2*sqrt((r1 - ra)^2 + (za - zj)^2)*((4*r1*ra/((r1 - ra)^2 + (za - zj)^2) + 1)*elliptic_kc(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2)) - elliptic_ec(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2)))*cI0*(2*(r1 - ra)*r1*ra/((r1 - ra)^2 + (za - zj)^2)^2 + r1/((r1 - ra)^2 + (za - zj)^2))/(r1*(4*r1*ra/((r1 - ra)^2 + (za - zj)^2) + 1)*ra) - 2*sqrt((r2 - ra)^2 + (za - zj)^2)*((4*r2*ra/((r2 - ra)^2 + (za - zj)^2) + 1)*elliptic_kc(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2)) - elliptic_ec(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2)))*cI0*(2*(r2 - ra)*r2*ra/((r2 - ra)^2 + (za - zj)^2)^2 + r2/((r2 - ra)^2 + (za - zj)^2))/(r2*(4*r2*ra/((r2 - ra)^2 + (za - zj)^2) + 1)*ra), zj, z1, z2)

# AS_diff_ra  =  - integrate(Iphi_js_rj.diff(ra).substitute(rj==r1), (zj, z1, z2), algorithm="giac") + integrate(Iphi_js_rj.diff(ra).substitute(rj==r2), (zj, z1, z2), algorithm="giac")
# AS_diff_ra  = -integrate(4*cI0*(r1 - ra)*elliptic_kc(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2))/((r1 - ra)^2 + (za - zj)^2)^(3/2) - 2*sqrt((r1 - ra)^2 + (za - zj)^2)*((4*r1*ra/((r1 - ra)^2 + (za - zj)^2) + 1)*elliptic_kc(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2)) - elliptic_ec(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2)))*cI0*(2*(r1 - ra)*r1*ra/((r1 - ra)^2 + (za - zj)^2)^2 + r1/((r1 - ra)^2 + (za - zj)^2))/(r1*(4*r1*ra/((r1 - ra)^2 + (za - zj)^2) + 1)*ra), zj, z1, z2) + integrate(4*cI0*(r2 - ra)*elliptic_kc(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2))/((r2 - ra)^2 + (za - zj)^2)^(3/2) - 2*sqrt((r2 - ra)^2 + (za - zj)^2)*((4*r2*ra/((r2 - ra)^2 + (za - zj)^2) + 1)*elliptic_kc(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2)) - elliptic_ec(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2)))*cI0*(2*(r2 - ra)*r2*ra/((r2 - ra)^2 + (za - zj)^2)^2 + r2/((r2 - ra)^2 + (za - zj)^2))/(r2*(4*r2*ra/((r2 - ra)^2 + (za - zj)^2) + 1)*ra), zj, z1, z2)

# Av_diff_ra = integrate(integrate(Iphi_jv_rj.diff(ra), (zj, z1, z2), algorithm="giac"), (rj, r1, r2), algorithm="giac")
# print "Av_diff_ra  =", Av_diff_ra
Av_diff_ra  = 0

#H_phi = AT_diff_za - (AS_diff_ra - Av_diff_ra)
#print "H_phi  =", H_phi
# H_phi  = integrate(4*cI0*(r1 - ra)*elliptic_kc(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2))/((r1 - ra)^2 + (za - zj)^2)^(3/2) - 2*sqrt((r1 - ra)^2 + (za - zj)^2)*((4*r1*ra/((r1 - ra)^2 + (za - zj)^2) + 1)*elliptic_kc(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2)) - elliptic_ec(-4*r1*ra/((r1 - ra)^2 + (za - zj)^2)))*cI0*(2*(r1 - ra)*r1*ra/((r1 - ra)^2 + (za - zj)^2)^2 + r1/((r1 - ra)^2 + (za - zj)^2))/(r1*(4*r1*ra/((r1 - ra)^2 + (za - zj)^2) + 1)*ra), zj, z1, z2) - integrate(4*cI0*(r2 - ra)*elliptic_kc(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2))/((r2 - ra)^2 + (za - zj)^2)^(3/2) - 2*sqrt((r2 - ra)^2 + (za - zj)^2)*((4*r2*ra/((r2 - ra)^2 + (za - zj)^2) + 1)*elliptic_kc(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2)) - elliptic_ec(-4*r2*ra/((r2 - ra)^2 + (za - zj)^2)))*cI0*(2*(r2 - ra)*r2*ra/((r2 - ra)^2 + (za - zj)^2)^2 + r2/((r2 - ra)^2 + (za - zj)^2))/(r2*(4*r2*ra/((r2 - ra)^2 + (za - zj)^2) + 1)*ra), zj, z1, z2) + integrate(4*cI0*(z1 - za)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2))/((ra - rj)^2 + (z1 - za)^2)^(3/2) - 4*((4*ra*rj/((ra - rj)^2 + (z1 - za)^2) + 1)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2)) - elliptic_ec(-4*ra*rj/((ra - rj)^2 + (z1 - za)^2)))*cI0*(z1 - za)/(((ra - rj)^2 + (z1 - za)^2)^(3/2)*(4*ra*rj/((ra - rj)^2 + (z1 - za)^2) + 1)), rj, r1, r2) - integrate(4*cI0*(z2 - za)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2))/((ra - rj)^2 + (z2 - za)^2)^(3/2) - 4*((4*ra*rj/((ra - rj)^2 + (z2 - za)^2) + 1)*elliptic_kc(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2)) - elliptic_ec(-4*ra*rj/((ra - rj)^2 + (z2 - za)^2)))*cI0*(z2 - za)/(((ra - rj)^2 + (z2 - za)^2)^(3/2)*(4*ra*rj/((ra - rj)^2 + (z2 - za)^2) + 1)), rj, r1, r2)

Z1 = -1
Z2 = 1
R1 = 1
R2 = 2
At_diff_za_substituted = At_diff_za.substitute(cI0==1, z1==Z1, z2==Z2)
print "At_diff_za_substituted =", At_diff_za_substituted

As_diff_ra_substituted = As_diff_ra.substitute(cI0==1, r1==R1, r2==R2)
print "As_diff_ra_substituted =", As_diff_ra_substituted

Za = -2
Ra = +1.5
At_diff_za_substituted2 = At_diff_za_substituted.substitute(za==Za, ra==Ra)
print "At_diff_za_substituted2 =", At_diff_za_substituted2
As_diff_ra_substituted2 = As_diff_ra_substituted.substitute(za==Za, ra==Ra)
print "As_diff_ra_substituted2 =", As_diff_ra_substituted2

At_diff_za_num_int = At_diff_za_substituted2.nintegral(rj, R1, R2)
print "At_diff_za_num_int  =", At_diff_za_num_int

As_v_diff_ra_num_int = (As_diff_ra_substituted2 - Av_diff_ra  ).nintegral(zj, Z1, Z2)
print "As_v_diff_ra_num_int  =", As_v_diff_ra_num_int

H_phi = At_diff_za_num_int[0] - As_v_diff_ra_num_int[0]
print "H_phi  =", H_phi

#At_diff_za_num_int  = (0.892985470996911, 9.91413030556511e-15, 21, 0)
#As_v_diff_ra_num_int  = (0.7922729550346237, 1.915043453188808e-14, 21, 0)
#H_phi  = 0.100712515962

import numpy as np

step_Ra = 0.5
min_Ra = -2.5
max_Ra = 2.5

ra_range = [min_Ra, max_Ra]
n_ra = int((max_Ra-min_Ra)/step_Ra) + 1

step_Za = 0.5
min_Za = -5
max_Za = 5

za_range = [min_Za, max_Za]
n_za = int((max_Za-min_Za)/step_Za) + 1



plot_data = []
m = []
for Ra in np.arange(min_Ra, max_Ra, step_Ra):
    h = []
    for Za in np.arange(min_Za, max_Za, step_Za):
        try:
            At_diff_za_substituted2 = At_diff_za_substituted.substitute(za==Za, ra==Ra)
            #print "At_diff_za_substituted2 =", At_diff_za_substituted2
            As_diff_ra_substituted2 = As_diff_ra_substituted.substitute(za==Za, ra==Ra)
            #print "As_diff_ra_substituted2 =", As_diff_ra_substituted2

            At_diff_za_num_int = At_diff_za_substituted2.nintegral(rj, R1, R2)
            #print "At_diff_za_num_int  =", At_diff_za_num_int

            As_v_diff_ra_num_int = (As_diff_ra_substituted2 - Av_diff_ra  ).nintegral(zj, Z1, Z2)
            #print "As_v_diff_ra_num_int  =", As_v_diff_ra_num_int

            H_phi = At_diff_za_num_int[0] - As_v_diff_ra_num_int[0]
            print "Ra  =", Ra, "Za  =", Za, "H_phi  =", H_phi

            plot_data.append((Za, Ra, H_phi))
            h.append(H_phi)
        except:
            pass
    m.append(h)

results_folder = "../articles/sagemath/field_of_deyna_cylinder/results/"
#results_folder = "./results/"

print plot_data
print m
_data = np.array(m, dtype = np.float64, copy = False, subok=True, ndmin = 0)

from sage.plot.contour_plot import ContourPlot
p = ContourPlot(m, za_range, ra_range, options=dict(fill=False, contours=[-1,0,1]))
g = Graphics()
g.add_primitive(p)

pname = results_folder + "H_phi_contour" + ".png"
print pname
g.save(pname)


g = list_plot3d(plot_data, frame_aspect_ratio=[1, 1, 1/3])


pname = results_folder + "H_phi" + ".png"
print pname
g.save(pname)
















