import sys
reload(sys)
sys.setdefaultencoding('utf8')

theta_0 = var("theta_0")
v = var("v")
c = var("c")
dt_zap = var("dt_zap")

x_1 = solve(x^2  +(x*cot(theta_0) + v*dt_zap)^2 - (c*dt_zap)^2, x)
print x_1

# sage: [
# x == -(dt_zap*v*cot(theta_0) + sqrt(c^2*cot(theta_0)^2 + c^2 - v^2)*dt_zap)/(cot(theta_0)^2 + 1),
# x == -(dt_zap*v*cot(theta_0) - sqrt(c^2*cot(theta_0)^2 + c^2 - v^2)*dt_zap)/(cot(theta_0)^2 + 1)
# ]

# sage: [
# x == -(dt_zap*v*cot(theta_0) + sqrt(c^2*cot(theta_0)^2 + c^2 - v^2)*dt_zap)/(cot(theta_0)^2 + 1),
# x == -(dt_zap*v*cot(theta_0) - sqrt(c^2*cot(theta_0)^2 + c^2 - v^2)*dt_zap)/(cot(theta_0)^2 + 1)
# ]

x1 = -(dt_zap*v*cot(theta_0) + sqrt(c^2*cot(theta_0)^2 + c^2 - v^2)*dt_zap)/(cot(theta_0)^2 + 1)

n = simplify(x1/(v*dt_zap))
print n

# -(dt_zap*v*cot(theta_0) + sqrt(c^2*cot(theta_0)^2 + c^2 - v^2)*dt_zap)/((cot(theta_0)^2 + 1)*dt_zap*v)
