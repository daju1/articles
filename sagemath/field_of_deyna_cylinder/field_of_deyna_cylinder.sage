import sys
reload(sys)
sys.setdefaultencoding('utf8')


z__j = var("z__j")
z__a = var("z__a")
rho__j = var("rho__j")
rho__a = var("rho__a")
varphi__j = var("varphi__j")

Ivarphi__j = integral(1/sqrt(rho__j^2+rho__a^2+(z__j-z__a)^2-2*rho__j*rho__a*cos(varphi__j)), varphi__j)
