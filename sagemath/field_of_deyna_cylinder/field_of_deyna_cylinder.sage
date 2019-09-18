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

integral(1/sqrt(rj^2+ra^2+(zj-za)^2-2*rj*ra*cos(phi)), (phi, 0, 2*pi))
# Ivarphi__j = integral(1/sqrt(rho__j^2+rho__a^2+(z__j-z__a)^2-2*rho__j*rho__a*cos(varphi__j)), (varphi__j,0,2*pi))
# RuntimeError: ECL says: Error executing code in Maxima:

rja2 = (rj-ra)^2+(zj-za)^2