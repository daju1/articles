#!/usr/bin/env python
# coding: utf-8

# О возникновении силы тяги в системе вращающихся зарядов

# А.Ю.Дроздов

# Пусть две массы (заряда или магнитных диполя) синхронно вращаются по окружностям с несовпадающим центром вращения в одном и том же направлении


# рассмотрим модель вращающихся зарядов и 
# рассчитаем поля движущихся по заданным траекториям зарядов 
# в соответствии с потенциалами Лиенара-Вихерта с учётом запаздывания

R = var("R")
S = var("S")
omega = var("omega")
t = var("t")
c = var("c")


# http://www.sciteclibrary.ru/cgi-bin/yabb2/YaBB.pl?num=1528093569/326#326
# 
# вот тот же вывод в более компактном виде. а если что-то непонятно, этот момент можно посмотреть в развернутой версии выше
# 
# $\vec{E} = -\nabla\varphi - \frac{1}{c}\frac{\partial}{\partial t}\vec{A} = \frac{q}{k^2}(\nabla k+\frac{\vec{v}}{c^2}\frac{\partial}{\partial t}k - \frac{\vec{a}r}{c^2})$, где $k = r - \frac{\vec{v}\vec{r}}{c}$
# 
# $dk = dr - \vec{r}d\vec{v}/c - \vec{v}d\vec{r}/c = (-c - \frac{\vec{r}\vec{a}}{c} + \frac{v^2}{c}) dt_1 + cdt_2 - \frac{\vec{v}}{c}d\vec{r_2}$
# 
# $\vec{E} = \frac{q}{k^2}(\nabla k+\frac{\vec{v}}{c^2}\frac{\partial}{\partial t}k - \frac{\vec{a}r}{c^2})$
# 
# $\vec{E} = \frac{q}{k^2}((-\frac{\vec{r}}{c k})(-c - \frac{\vec{r}\vec{a}}{c} + \frac{v^2}{c})-\frac{\vec{v}}{c}+\frac{\vec{v}}{c^2}(\frac{r}{k}(-c - \frac{\vec{r}\vec{a}}{c} + \frac{v^2}{c}) + c) - \frac{\vec{a}r}{c^2})$
# 
# $\vec{E} = \frac{q}{k^3}((\vec{r}-\frac{r}{c}\vec{v})(1 + \frac{\vec{r}\vec{a}}{c^2} - \frac{v^2}{c^2}) - \vec{a}\frac{kr}{c^2})$
# 
# и магнитное поле аналогичным образом:
# 
# $\frac{\partial}{\partial x}\vec{A} = \frac{\partial}{\partial x}\frac{q\vec{v}}{k c}$
# 
# $\frac{\partial}{\partial x}\vec{A} = \frac{q}{k c}\frac{\partial}{\partial x} \vec{v} - \frac{q\vec{v}}{k^2 c}\frac{\partial}{\partial x} k$
# 
# $\frac{\partial}{\partial x}\vec{A} = \frac{q}{k c}(-\frac{\vec{a}r_x}{k c})- \frac{q\vec{v}}{k^2 c}(\frac{r_x}{k}(1 +\frac{\vec{r}\vec{a}}{c^2} - \frac{v^2}{c^2}) - \frac{v_x}{c})$
# 
# $\vec{B} = \nabla\times\vec{A} = -(\frac{q}{k^2 c^2})\vec{r}\times\vec{a} - \frac{q}{k^3 c}(1 + \frac{\vec{r}\vec{a}}{c^2} - \frac{v^2}{c^2})\vec{r}\times\vec{v}$
# 
# $\vec{B} = -\frac{q}{k^3}(\vec{r}\times\frac{\vec{v}}{c}(1 + \frac{\vec{r}\vec{a}}{c^2} - \frac{v^2}{c^2}) + \vec{r}\times\frac{\vec{a}k}{c^2})$
# 
# $\vec{B} = \frac{\vec{r}\times\vec{E}}{r}$
# 

# http://www.sciteclibrary.ru/cgi-bin/yabb2/YaBB.pl?num=1528093569/413#413
# 
# вот тот же вывод в более компактном виде. а если что-то непонятно, этот момент можно посмотреть в развернутой версии выше
# 
# $\vec{E} = \frac{q}{k^3}((\vec{r}-\frac{r}{c}\vec{v})(1 + \frac{\vec{r}\vec{a}}{c^2} - \frac{v^2}{c^2})) - \vec{a}\frac{kr}{c^2})$
# 
# 
# Что-то в этой формуле не так: открывающихся скобок три, а закрывающихся четыре...
# 
# 
# 
# тут вроде и из размерности и из предыдущей в выводе видно где скобки:
# 
# 
# 
# $\vec{E} = \frac{q}{k^2}(\nabla k+\frac{\vec{v}}{c^2}\frac{\partial}{\partial t}k - \frac{\vec{a}r}{c^2})$
# 
# $\vec{E} = \frac{q}{k^2}((-\frac{\vec{r}}{c k})(-c - \frac{\vec{r}\vec{a}}{c} + \frac{v^2}{c})-\frac{\vec{v}}{c}+\frac{\vec{v}}{c^2}(\frac{r}{k}(-c - \frac{\vec{r}\vec{a}}{c} + \frac{v^2}{c}) + c) - \frac{\vec{a}r}{c^2})$
# 
# $\vec{E} = \frac{q}{k^3}((\vec{r}-\frac{r}{c}\vec{v})(1 + \frac{\vec{r}\vec{a}}{c^2} - \frac{v^2}{c^2}) - \vec{a}\frac{kr}{c^2})$
# 
# 

# In[6]:


# расчет итерациями запаздывающего момента
Epsilon = 10^-14 # погрешность
no_retardition_test = False
def tlag(x, y, z, t, s):
    if no_retardition_test:
        return (t, s(t)[0], s(t)[1], s(t)[2])
    # print(s)
    t1 = t
    t2 = t - 2 * Epsilon
    #print(s(t1))
    while (abs(t1 - t2) > Epsilon):
        t1 = t2
        t2 = t - sqrt((x - s(t1)[0])^2 + (y - s(t1)[1])^2 + (z - s(t1)[2])^2) / c
        #print(t2)
    return (t2, s(t2)[0], s(t2)[1], s(t2)[2])


# Здесь нужно отметить во избежание путаницы, что радиус Лиенара Вихерта $k$ в формуле для электрического и магнитного полей и $k$ в программе отличаются друг от друга тем, что в программе $k$ нормирован на единицу, а в формуле нет. При переходе от формул к програмным кодам по сути произведено преобразование $k \rightarrow k\cdot r$

# In[7]:


def calc_k(x, y, z, t, s, v, t2):
    if no_retardition_test:
        r = sqrt((x - s(t)[0])^2 + (y - s(t)[1])^2 + (z - s(t)[2])^2)
    else:
        r = c * (t - t2)
    nx = (x - s(t2)[0])/r
    ny = (y - s(t2)[1])/r
    nz = (z - s(t2)[2])/r
    if no_retardition_test:
        k = 1.0
    else:
        k = 1.0 - (nx*v(t2)[0] + ny * v(t2)[1] + nz * v(t2)[2]) / c
    return k, r, nx, ny, nz


# In[8]:


# отношение радиуса Лиенара Вихерта к радиусу
def klw(x, y, z, t, s, v, tlag=tlag):
    # = Module[{t2, r, nx, k},
    (t2, x_, y_, z_) = tlag(x, y, z, t, s) # расчет итерациями запаздывающего момента
    (k, r, nx, ny, nz) = calc_k(x, y, z, t, s, v, t2)
    return k


# In[9]:


# Радиус Лиенара Вихерта
def Rlw(x, y, z, t, s, v, tlag=tlag):
    # Module[{t2, r, nx, k},
    (t2, x_, y_, z_) = tlag(x, y, z, t, s) # расчет итерациями запаздывающего момента
    (k, r, nx, ny, nz) = calc_k(x, y, z, t, s, v, t2)
    return k * r


# http://www.sciteclibrary.ru/cgi-bin/yabb2/YaBB.pl?num=1528093569/339#339
# 
# чему равен
# $\nabla k$, где $k = r - \frac{\vec{v}\vec{r}}{c}$
# а также чему равен 
# $\nabla \frac{1}{k}$
# 
# 
# $\nabla k = (\frac{\vec{r}}{k})(1 + \frac{\vec{r}\vec{a}}{c^2}-\frac{v^2}{c^2})-\frac{\vec{v}}{c}$
# 
# $\nabla\frac{1}{k} = (\frac{d}{dk}\frac{1}{k})(\nabla k) = -\frac{1}{k^2}\nabla k$

# In[10]:


# phi - скалярный потенциал Лиенара Вихерта
def phi(x, y, z, t, s, v, tlag=tlag):
    # phi[x_, y_, t_] := Module[{t2, r, nx, k},
    (t2, x_, y_, z_) = tlag(x, y, z, t, s) # расчет итерациями запаздывающего момента
    (k, r, nx, ny, nz) = calc_k(x, y, z, t, s, v, t2)
    return 1.0/(k*r)


# $\vec{E} = \frac{q}{k^2}(\nabla k+\frac{\vec{v}}{c^2}\frac{\partial}{\partial t}k - \frac{\vec{a}r}{c^2})$
# 
# $\vec{E} = \frac{q}{k^2}((-\frac{\vec{r}}{c k})(-c - \frac{\vec{r}\vec{a}}{c} + \frac{v^2}{c})-\frac{\vec{v}}{c}+\frac{\vec{v}}{c^2}(\frac{r}{k}(-c - \frac{\vec{r}\vec{a}}{c} + \frac{v^2}{c}) + c) - \frac{\vec{a}r}{c^2})$
# 

# In[11]:


# эл. поле - grad phi
def electrgradphi(x, y, z, t, s, v, w, tlag=tlag):
    # electrgradphi[x_, y_, t_] := Module[{t2, r, nx, ny, n, k},
    (t2, x_, y_, z_) = tlag(x, y, z, t, s) # расчет итерациями запаздывающего момента
    (k, r, nx, ny, nz) = calc_k(x, y, z, t, s, v, t2)
    v2_c2 = (v(t2)[0]^2 + v(t2)[1]^2 + v(t2)[2]^2) / c^2
    ra_c2 = r * (nx*w(t2)[0] + ny*w(t2)[1] + nz*w(t2)[2]) / c^2
    #print ("v2_c2", v2_c2)
    #print ("ra_c2", ra_c2)
    #print ("v(t2)", v(t2))
    E1_x = (1.0/(r^2*k^2))*((1 - v2_c2 + ra_c2) * nx/k - v(t2)[0]/c)
    E1_y = (1.0/(r^2*k^2))*((1 - v2_c2 + ra_c2) * ny/k - v(t2)[1]/c)
    E1_z = (1.0/(r^2*k^2))*((1 - v2_c2 + ra_c2) * nz/k - v(t2)[2]/c)
    return (E1_x, E1_y, E1_z)


# In[12]:


# эл. поле dA_dt
def electrdAdt(x, y, z, t, s, v, w, tlag=tlag):
    # electrdAdt[x_, y_, t_] := Module[{t2, r, nx, ny, n, k},
    (t2, x_, y_, z_) = tlag(x, y, z, t, s) # расчет итерациями запаздывающего момента
    (k, r, nx, ny, nz) = calc_k(x, y, z, t, s, v, t2)
    # n = (nx,ny,nz)
    v2_c2 = (v(t2)[0]^2 + v(t2)[1]^2 + v(t2)[2]^2) / c^2
    ra_c2 = r * (nx*w(t2)[0] + ny*w(t2)[1] + nz*w(t2)[2]) / c^2
    E2_x = (1.0/(r^2*k^2))*(v(t2)[0]/c * (1/k * (v2_c2 - ra_c2 - 1.0) + 1.0) - r*w(t2)[0]/c^2)
    E2_y = (1.0/(r^2*k^2))*(v(t2)[1]/c * (1/k * (v2_c2 - ra_c2 - 1.0) + 1.0) - r*w(t2)[1]/c^2)
    E2_z = (1.0/(r^2*k^2))*(v(t2)[2]/c * (1/k * (v2_c2 - ra_c2 - 1.0) + 1.0) - r*w(t2)[2]/c^2)
    return (E2_x, E2_y, E2_z)


# $\vec{E} = \frac{q}{k^3}\left(\left(\vec{r}-\frac{r}{c}\vec{v}\right)\left(1 + \frac{\vec{r}\vec{a}}{c^2} - \frac{v^2}{c^2}\right) - \vec{a}\frac{kr}{c^2}\right)$
# 
# При переходе от формул к програмным кодам произведено преобразование
# 
# $k \rightarrow k\cdot r$
# 
# $\vec{r} \rightarrow r\cdot \vec{n}$
# 
# $\vec{E} = \frac{q}{k^3}\left(\left(\vec{n}-\frac{\vec{v}}{c}\right)\left(1 + \frac{\vec{r}\vec{a}}{c^2} - \frac{v^2}{c^2}\right) \frac{1}{r^2} - \frac{k}{r}\frac{\vec{a}}{c^2}\right)$

# In[13]:


# эл. поле
def electr(x, y, z, t, s, v, w, tlag=tlag):
    # electr[x_, y_, t_] := Module[{t2, r, nx, ny, n, k},
    (t2, x_, y_, z_) = tlag(x, y, z, t, s) # расчет итерациями запаздывающего момента
    (k, r, nx, ny, nz) = calc_k(x, y, z, t, s, v, t2)
    n = (nx,ny,nz)
    v2_c2 = (v(t2)[0]^2 + v(t2)[1]^2 + v(t2)[2]^2) / c^2
    ra_c2 = r * (nx*w(t2)[0] + ny*w(t2)[1] + nz*w(t2)[2]) / c^2
    E_x = (1.0/k^3)*((1 - v2_c2 + ra_c2)*(nx - v(t2)[0]/c)/r^2 - (k/r)*w(t2)[0]/c^2)
    E_y = (1.0/k^3)*((1 - v2_c2 + ra_c2)*(ny - v(t2)[1]/c)/r^2 - (k/r)*w(t2)[1]/c^2)
    E_z = (1.0/k^3)*((1 - v2_c2 + ra_c2)*(nz - v(t2)[2]/c)/r^2 - (k/r)*w(t2)[2]/c^2)
    return (E_x, E_y, E_z, x_, y_, z_)


# $\vec{B} = \frac{\vec{r}\times\vec{E}}{r}$
# 
# $\vec{B} = -\frac{q}{k^3}\left(\vec{r}\times\frac{\vec{v}}{c}\left(1 + \frac{\vec{r}\vec{a}}{c^2} - \frac{v^2}{c^2}\right) + \vec{r}\times\frac{\vec{a}k}{c^2}\right)$
# 
# При переходе от формул к програмным кодам произведено преобразование
# 
# $k \rightarrow k\cdot r$
# 
# $\vec{r} \rightarrow r\cdot \vec{n}$
# 
# $\vec{B} = -\frac{q}{k^3}\left(
# \left(\vec{n}\times{\vec{v}}\right)
# \left(1 + \frac{\vec{r}\vec{a}}{c^2} - \frac{v^2}{c^2}\right)  \frac{1}{c\,r^2} + \left(\vec{n}\times\vec{a}\right)
# \frac{k}{r\,c^2}\right)$

# $B_x = -\frac{q}{k^3}\left(
# \left(   n_y\,v_z -  n_z\,v_y \right)
# \left(1 + \frac{\vec{r}\vec{a}}{c^2} - \frac{v^2}{c^2}\right)  \frac{1}{c\,r^2} +
# \left(  n_y\,a_z -  n_z\,a_y  \right)
# \frac{k}{r\,c^2}\right)$
# 
# $B_y = -\frac{q}{k^3}\left(
# \left(   n_z\,v_x -  n_x\,v_z \right)
# \left(1 + \frac{\vec{r}\vec{a}}{c^2} - \frac{v^2}{c^2}\right)  \frac{1}{c\,r^2} +
# \left(  n_z\,a_x -  n_x\,a_z  \right)
# \frac{k}{r\,c^2}\right)$
# 
# $B_z = -\frac{q}{k^3}\left(
# \left(   n_x\,v_y -  n_y\,v_x \right)
# \left(1 + \frac{\vec{r}\vec{a}}{c^2} - \frac{v^2}{c^2}\right)  \frac{1}{c\,r^2} +
# \left(  n_x\,a_y -  n_y\,a_x  \right)
# \frac{k}{r\,c^2}\right)$

# In[14]:


# магнитное поле
def magnet(x, y, z, t, s, v, w, tlag=tlag):
    # electr[x_, y_, t_] := Module[{t2, r, nx, ny, n, k},
    (t2, x_, y_, z_) = tlag(x, y, z, t, s) # расчет итерациями запаздывающего момента
    (k, r, nx, ny, nz) = calc_k(x, y, z, t, s, v, t2)
    n = (nx,ny,nz)
    v2_c2 = (v(t2)[0]^2 + v(t2)[1]^2 + v(t2)[2]^2) / c^2
    ra_c2 = r * (nx*w(t2)[0] + ny*w(t2)[1] + nz*w(t2)[2]) / c^2
    B_x = -(1.0/k^3)*((1 - v2_c2 + ra_c2)*(ny*v(t2)[2] - nz*v(t2)[1])/r^2/c + (ny*w(t2)[2] - nz*w(t2)[1])*(k/r)/c^2)
    B_y = -(1.0/k^3)*((1 - v2_c2 + ra_c2)*(nz*v(t2)[0] - nx*v(t2)[2])/r^2/c + (nz*w(t2)[0] - nx*w(t2)[2])*(k/r)/c^2)
    B_z = -(1.0/k^3)*((1 - v2_c2 + ra_c2)*(nx*v(t2)[1] - ny*v(t2)[0])/r^2/c + (nx*w(t2)[1] - nx*w(t2)[0])*(k/r)/c^2)
    return (B_x, B_y, B_z, x_, y_, z_)


# In[15]:


# Drobyshev's test
# Как запаздывающий Лиенар-Вихерт становится "незапаздывающим". Визуализация
# http://www.sciteclibrary.ru/cgi-bin/yabb2/YaBB.pl?num=1528093569/0

# c = 1
# vk = 0.84 # finish velocity
# a = 0.3   # acseleration
# t0 = vk/a # time of acseleration

# s = lambda t : (0 if t < 0 else (a * t^2 / 2 if t < t0 else vk * t - a * t0^2 / 2), 0, 0)
# v = lambda t : (0 if t < 0 else (a * t       if t < t0 else vk)                   , 0, 0)
# w = lambda t : (0 if t < 0 else (a           if t < t0 else 0)                    , 0, 0)


# In[19]:


#sgs 

omega_d = 2*pi.n()*150000/60
omega_d *= 10^5
c = 299792458 * 100
R_r = 10
R_l = 9


# In[17]:


print("omega_d =", omega_d)
print(omega_d * R_r / c)
print(omega_d * R_l / c)

# In[18]:


S = 0.99

# centers of circles
cr = (0, 0, 0)
cl = (S/2, 0, 0)


# In[19]:


# координата от времени
omega = var("omega")
#c = var("c")
# current positions of rotated masses
pr = lambda t : (cr[0] + R_r*cos(omega * t), cr[1] + R_r*sin(omega * t), 0)
pl = lambda t : (cl[0] + R_l*cos(omega * t), cl[1] + R_l*sin(omega * t), 0)

print(pr(t))
print(pl(t))


# In[20]:


# скорость от времени
exec(preparse("vr = lambda t : " + str((diff(pr(t)[0],t), diff(pr(t)[1],t), diff(pr(t)[2],t)))))
exec(preparse("vl = lambda t : " + str((diff(pl(t)[0],t), diff(pl(t)[1],t), diff(pl(t)[2],t)))))

print(vr(t))
print(vl(t))


# In[21]:


# ускорение
exec(preparse("ar = lambda t : " + str((diff(vr(t)[0],t), diff(vr(t)[1],t), diff(vr(t)[2],t)))))
exec(preparse("al = lambda t : " + str((diff(vl(t)[0],t), diff(vl(t)[1],t), diff(vl(t)[2],t)))))

print(ar(t))
print(al(t))


# In[22]:

def find_newton_root(f,x,xn, max_steps=20, debug = False):
    if debug:
        print("find_newton_root f =", f)
    df = f.diff(x)
    if debug:
        print ("df", df)
        print ("f/df", f/df)
    step = 1.0
    #NewtonIt = lambda x_, step : x_-step*(f/df).subs(x == x_)
    #xn=(a+b)/2;                      # initial guess
    
    def NewtonIt(_x, step):
        if debug:
            print("_x", _x)
        delta = (f/df).subs(x == _x)
        if debug:
            print("delta", delta)
            print("step", step)
        step_delta = step*delta
        if debug:
            print("step_delta", step_delta)

        res = _x-step_delta
        if debug:
            print ("_x-step_delta", _x-step_delta)
            print ("res", res)
        return res

    if debug:
        print ("xn", xn)
    for i in range(max_steps):
        #xn=N(NewtonIt(xn, step), digits=32)\n",
        if float != type(xn):
            xn = xn.n()
        xn=NewtonIt(xn, step)
        if debug:
            print ("xn", xn)

        if debug:
            f_n = f.subs(x == xn)
            print ("f_n", f_n)
        #step *= 0.999

    return xn


# In[23]:


logging = False
tol=0.000000000001
def find_root_recursive(func,a,b,tol=0.000000000001):
    try:
        free_variable = func.variables()[0]
    except:
        free_variable = x
    #print("free_variable", free_variable)
    #print("func", func(1))
    L = []
    if b - a < tol:
        return L
    try:
        #print ("a = ", a, "b = ", b)
        #print (func(free_variable=a), func(free_variable=b))
        x0 = find_root(func,a,b)
        #exec("print(func(" + preparse(str(free_variable)) + "=x0))")
        #print ("x0 =", x0, func(free_variable=x0))
        L.append(x0)
        L += find_root_recursive(func,a,x0-tol,tol)       
        L += find_root_recursive(func,x0+tol,b,tol)       
    except Exception as ex:
        if 'f appears to have no zero on the interval' != str(ex):
            if logging:
                print(str(ex))
                print ("a = ", a, "b = ", b)
                print ("func = ", func)
                print ("func(a) =", func.subs(free_variable==a))
                print ("func(b) =", func.subs(free_variable==b))
                #exec("print(func(" + preparse(str(free_variable)) + "=a))")
                #exec("print(func(" + preparse(str(free_variable)) + "=b))")
        pass
    L.sort()
    return L


# In[24]:


# расчет итерациями запаздывающего момента

no_retardition_test = False
def tlag_symbolic(x, y, z, t, s):
    if no_retardition_test:
        return (t, s(t)[0], s(t)[1], s(t)[2])
    tau = var('tau')
    # eq = c^2 * (t-tau)^2 == (x - s(tau)[0])^2 + (y - s(tau)[1])^2 + (z - s(tau)[2])^2
    eq_rhs = (x - s(tau)[0])^2 + (y - s(tau)[1])^2 + (z - s(tau)[2])^2
    eq_rhs = eq_rhs.expand().subs(cos(omega*tau)^2 == 1 - sin(omega*tau)^2)
    #print("eq_rhs=", eq_rhs)
    eq = c^2 * (t-tau)^2 == eq_rhs
    solved = False
    if not eq_rhs.has(tau):
        t2_symbolic = solve(eq, tau)
        if logging:
            print("Logging: tlag")
            print("s =", s(tau))
            print("eq =", eq)
            print("sol =", t2_symbolic)

        for t2_symb in t2_symbolic:
            t2_i = t2_symb.rhs()
            t2_i = t2_i.subs(cos(omega*tau)^2 + sin(omega*tau)^2 == 1)
            if logging:
                print ("t2_i=", t2_i)
            if not t2_i.has(tau):
                if t2_i < t:
                    t2 = t2_i
                    solved = True

    if False == solved:
        if logging:
            print("Warning: tlag")
            print("s =", s(tau))
            print("eq =", eq)
        eq_subs = eq.subs(omega == omega_d)
        if logging:
            print("eq_subs =", eq_subs)
        delta_tau = (2*pi / omega_d).n()
        if logging:
            print("delta_tau =", delta_tau)
        a = t-delta_tau
        b = t
        if logging:
            print("a,b =", a, ",", b)
        eq_tau = (eq_subs.lhs()-eq_subs.rhs())#.subs(tau == var("x"))
        if logging:
            print("eq_tau =", eq_tau)
        t2_roots = find_root_recursive(eq_tau, a, b)
        if logging:
            print("t2_roots =", t2_roots)
        if len(t2_roots) > 0:
            t2_root = t2_roots[len(t2_roots)-1]
            t2 = find_newton_root(f = eq_tau, x = tau, xn = t2_root)
        else:
            if logging:
                #plot(eq_tau, var("x"), a, b).show()
                plot(eq_tau, tau, a, b).show()
            eq_tau_a = eq_tau.subs(tau == a)
            eq_tau_b = eq_tau.subs(tau == b)

            n_try = 0
            while (eq_tau_a * eq_tau_b > 0 and n_try < 10):
                a -= delta_tau
                eq_tau_a = eq_tau.subs(tau == a)
                n_try += 1

            if (eq_tau_a * eq_tau_b < 0):
                t2 = find_newton_root(f = eq_tau, x = tau, xn = (a+b)/2, max_steps=50)
            else:
                print("eq_tau_a", eq_tau_a)
                print("eq_tau_b", eq_tau_b)
                plot(eq_tau, tau, a, b).show()
                raise Exception("t2 roots are empty for eq='" + str(eq) + "' eq_subs='" + str(eq_subs) + "'" + "in the interval [" + str(a) + "," + str(b) + "]")
            
    return (t2, s(t2)[0], s(t2)[1], s(t2)[2])


# In[ ]:


# In[30]:


(E1_x, E1_y, E1_z) = electrgradphi(0, 0, 0, 1, pr, vr, ar)
print(E1_x.subs(omega == omega_d))
print(E1_y.subs(omega == omega_d))
print(E1_z.subs(omega == omega_d))


# In[31]:


(E2_x, E2_y, E2_z) = electrdAdt(0, 0, 0, 1, pr, vr, ar)
print(E2_x.subs(omega == omega_d))
print(E2_y.subs(omega == omega_d))
print(E2_z.subs(omega == omega_d))


# In[32]:


(E1_x, E1_y, E1_z)          = electrgradphi(0, 0, 0, 1, pr, vr, ar)
(E2_x, E2_y, E2_z)          = electrdAdt(0, 0, 0, 1, pr, vr, ar)
(E_x, E_y, E_z, x_, y_, z_) = electr(0, 0, 0, 1, pr, vr, ar)
(B_x, B_y, B_z, x_, y_, z_) = magnet(0, 0, 0, 1, pr, vr, ar)

print(E_x.subs(omega == omega_d), E1_x.subs(omega == omega_d)+E2_x.subs(omega == omega_d))
print(E_y.subs(omega == omega_d), E1_y.subs(omega == omega_d)+E2_y.subs(omega == omega_d))
print(E_z.subs(omega == omega_d), E1_z.subs(omega == omega_d)+E2_z.subs(omega == omega_d))
print (B_x.subs(omega == omega_d), (B_y).subs(omega == omega_d), (B_z).subs(omega == omega_d))

# In[33]:


def get_integrand_view(f):
    return f(x)

class my_dummy_integral:
    f = None
    a = None
    b = None
    def __init__(self, f, a, b):
        print ("my_dummy_integral ", f, a, b)
        self.f = f
        self.a = a
        self.b = b

def num_int(f, a, b):
    from scipy import integrate

    to_call_integration = True

    if type(f) is my_dummy_integral:
        to_call_integration = False

    import inspect
    stack = inspect.stack()
    for frame in stack:
        func_name = frame[3]
        # print ("func_name = ", func_name)
        if ('get_integrand_view' == func_name):
            to_call_integration = False
            break;

    if not to_call_integration:
        return my_dummy_integral(f,a,b)

    try:
        integral = integrate.quad(f, a, b)

        result = integral[0]
        return result

    except Exception as ex:

        print ("Exception ex = ", str(ex))
        print ("f = ", f)
        try:
            print ("integrand = ", get_integrand_view(f))
        except Exception as ex2:
            print ("Exception ex2 = ", ex2)
        print ("a = ", a)
        print ("b = ", b)

        raise ex

        integral = numerical_integral(f, a, b)

        print ("integral = ", integral)

        result = integral[0]
        print ("result = ", result)
        return result


# Расчёт углового усилия

# In[35]:


def calc_F_alpha(E, alpha_R, sign_L_sign_R):
    f_x = sign_L_sign_R * E[0]
    f_y = sign_L_sign_R * E[1]
    return f_y * cos(alpha_R) - f_x * sin(alpha_R)


# In[36]:


def cals_sum_F(n, S, R_l, R_r, alpha0_l = 0, alpha0_r = 0, omega_l = omega_d, omega_r = omega_d, to_animate = True):
    
    t = var("t")

    # centers of circles
    cr = (S/2, 0, 0)
    cl = (0, 0, 0)

    # current positions of rotated masses
    sign_r = []
    sign_l = []
    alpha_r = []
    alpha_l = []
    p_r = []
    p_l = []
    v_r = []
    v_l = []
    a_r = []
    a_l = []

    for i in range(n):
        sign_r += [lambda i=i : -((i%2)*2-1)]
        sign_l += [lambda i=i :   (i%2)*2-1]

        alpha_r += [lambda t, i=i : omega_r * t + i * 2*pi.n()/n + alpha0_r]
        alpha_l += [lambda t, i=i : omega_l * t + i * 2*pi.n()/n + alpha0_l]
        # координата от времени
        p_r += [lambda t, i=i : (cr[0] + R_r*cos(alpha_r[i](t)), cr[1] + R_r*sin(alpha_r[i](t)), 0)]
        p_l += [lambda t, i=i : (cl[0] + R_l*cos(alpha_l[i](t)), cl[1] + R_l*sin(alpha_l[i](t)), 0)]
        # скорость от времени
        exec(preparse("v_r += [lambda t, i=i: " + str((diff(p_r[i](t)[0],t), diff(p_r[i](t)[1],t), diff(p_r[i](t)[2],t))) + "]"))
        exec(preparse("v_l += [lambda t, i=i: " + str((diff(p_l[i](t)[0],t), diff(p_l[i](t)[1],t), diff(p_l[i](t)[2],t))) + "]"))
        # ускорение
        exec(preparse("a_r += [lambda t, i=i: " + str((diff(v_r[i](t)[0],t), diff(v_r[i](t)[1],t), diff(v_r[i](t)[2],t))) + "]"))
        exec(preparse("a_l += [lambda t, i=i: " + str((diff(v_l[i](t)[0],t), diff(v_l[i](t)[1],t), diff(v_l[i](t)[2],t))) + "]"))
        
        if False:
            print("alpha_r[", i, "]", alpha_r[i](t))
            print("alpha_l[", i, "]", alpha_l[i](t))

            print("p_r[", i, "]", p_r[i](t))
            print("p_l[", i, "]", p_l[i](t))

            print("v_r[", i, "]", v_r[i](t))
            print("v_l[", i, "]", v_l[i](t))

            print("a_r[", i, "]", a_r[i](t))
            print("a_l[", i, "]", a_l[i](t))


    Fx_l = []
    Fy_l = []
    Fx_r = []
    Fy_r = []
    
    F_alpha_l = []
    F_alpha_r = []
    
    omega_m = min(omega_l, omega_r)

    for i_a in range(n):
        for i_q in range(n):
            i_l = i_a 
            i_r = i_q
            
            sign_a = sign_l[i_a]()
            sign_q = sign_r[i_q]()
            
            sign_L_sign_R = sign_a*sign_q

            #XYZa = p_l[i_a](t_i)
            #XYZq = p_r[i_q](t_i)
            #Xa = XYZa[0]
            #Ya = XYZa[1]
            #Za = XYZa[2]
            #Xq = XYZq[0]
            #Yq = XYZq[1]
            #Zq = XYZq[2]

            #print(Xa,Ya,Za)
            #(E_x, E_y, E_z, x_, y_, z_) = electr(Xa, Ya, Za, t_i, p_r[i_q], v_r[i_q], a_r[i_q], tlag=tlag_symbolic)
            E1_x = num_int(lambda t_i : electrgradphi(p_l[i_a](t_i)[0], p_l[i_a](t_i)[1], p_l[i_a](t_i)[2], t_i, p_r[i_q], v_r[i_q], a_r[i_q], tlag=tlag_symbolic)[0], 0, 2*pi.n()/omega_m)
            E1_y = num_int(lambda t_i : electrgradphi(p_l[i_a](t_i)[0], p_l[i_a](t_i)[1], p_l[i_a](t_i)[2], t_i, p_r[i_q], v_r[i_q], a_r[i_q], tlag=tlag_symbolic)[1], 0, 2*pi.n()/omega_m)
            E2_x = num_int(lambda t_i : electrdAdt   (p_l[i_a](t_i)[0], p_l[i_a](t_i)[1], p_l[i_a](t_i)[2], t_i, p_r[i_q], v_r[i_q], a_r[i_q], tlag=tlag_symbolic)[0], 0, 2*pi.n()/omega_m)
            E2_y = num_int(lambda t_i : electrdAdt   (p_l[i_a](t_i)[0], p_l[i_a](t_i)[1], p_l[i_a](t_i)[2], t_i, p_r[i_q], v_r[i_q], a_r[i_q], tlag=tlag_symbolic)[1], 0, 2*pi.n()/omega_m)
            fx_l = (1.0*E1_x+E2_x)*sign_a*sign_q
            fy_l = (1.0*E1_y+E2_y)*sign_a*sign_q
            
            #E_x = num_int(lambda t_i : electr        (p_l[i_a](t_i)[0], p_l[i_a](t_i)[1], p_l[i_a](t_i)[2], t_i, p_r[i_q], v_r[i_q], a_r[i_q], tlag=tlag_symbolic)[0], 0, 2*pi.n()/omega_m)
            #E_y = num_int(lambda t_i : electr        (p_l[i_a](t_i)[0], p_l[i_a](t_i)[1], p_l[i_a](t_i)[2], t_i, p_r[i_q], v_r[i_q], a_r[i_q], tlag=tlag_symbolic)[1], 0, 2*pi.n()/omega_m)
            #print("fx_l,Ex", fx_l,E_x)
            #print("fy_l,Ey", fy_l,E_y)
            
            F_alpha_L = num_int(lambda t_i : calc_F_alpha(electr(p_l[i_a](t_i)[0], p_l[i_a](t_i)[1], p_l[i_a](t_i)[2], t_i, p_r[i_q], v_r[i_q], a_r[i_q], tlag=tlag_symbolic), alpha_l[i_l](t_i), sign_L_sign_R), 0, 2*pi.n()/omega_d)

            
            Fx_l += [fx_l]
            Fy_l += [fy_l]
            F_alpha_l += [F_alpha_L]
            #X2_q = x_
            #Y2_q = y_

            #(E_x, E_y, E_z, x_, y_, z_) = electr(Xq, Yq, Zq, t_i, p_l[i_a], v_l[i_a], a_l[i_a], tlag=tlag_symbolic)
            E1_x = num_int(lambda t_i : electrgradphi(p_r[i_q](t_i)[0], p_r[i_q](t_i)[1], p_r[i_q](t_i)[2], t_i, p_l[i_a], v_l[i_a], a_l[i_a], tlag=tlag_symbolic)[0], 0, 2*pi.n()/omega_m)
            E1_y = num_int(lambda t_i : electrgradphi(p_r[i_q](t_i)[0], p_r[i_q](t_i)[1], p_r[i_q](t_i)[2], t_i, p_l[i_a], v_l[i_a], a_l[i_a], tlag=tlag_symbolic)[1], 0, 2*pi.n()/omega_m)
            E2_x = num_int(lambda t_i : electrdAdt   (p_r[i_q](t_i)[0], p_r[i_q](t_i)[1], p_r[i_q](t_i)[2], t_i, p_l[i_a], v_l[i_a], a_l[i_a], tlag=tlag_symbolic)[0], 0, 2*pi.n()/omega_m)
            E2_y = num_int(lambda t_i : electrdAdt   (p_r[i_q](t_i)[0], p_r[i_q](t_i)[1], p_r[i_q](t_i)[2], t_i, p_l[i_a], v_l[i_a], a_l[i_a], tlag=tlag_symbolic)[1], 0, 2*pi.n()/omega_m)
            fx_r = (1.0*E1_x+E2_x)*sign_a*sign_q
            fy_r = (1.0*E1_y+E2_y)*sign_a*sign_q
            
            #E_x = num_int(lambda t_i : electr        (p_r[i_q](t_i)[0], p_r[i_q](t_i)[1], p_r[i_q](t_i)[2], t_i, p_l[i_a], v_l[i_a], a_l[i_a], tlag=tlag_symbolic)[0], 0, 2*pi.n()/omega_m)
            #E_y = num_int(lambda t_i : electr        (p_r[i_q](t_i)[0], p_r[i_q](t_i)[1], p_r[i_q](t_i)[2], t_i, p_l[i_a], v_l[i_a], a_l[i_a], tlag=tlag_symbolic)[1], 0, 2*pi.n()/omega_m)
            #print("fx_r,Ex", fx_r,E_x)
            #print("fy_r,Ey", fy_r,E_y)
            
            F_alpha_R = num_int(lambda t_i : calc_F_alpha(electr(p_r[i_q](t_i)[0], p_r[i_q](t_i)[1], p_r[i_q](t_i)[2], t_i, p_l[i_a], v_l[i_a], a_l[i_a], tlag=tlag_symbolic), alpha_r[i_r](t_i), sign_L_sign_R), 0, 2*pi.n()/omega_m)
            
            Fx_r += [fx_r]
            Fy_r += [fy_r]
            F_alpha_r += [F_alpha_R]
            #X2_a = x_
            #Y2_a = y_

            #print("fx_l=%f fy_l=%f fx_r=%f fy_r=%f" % (fx_l, fy_l, fx_r, fy_r))

    Fxt_l = []
    Fyt_l = []
    Fxt_r = []
    Fyt_r = []
    frames = []
    if to_animate:
        for t_i in tlist:
            if to_animate:
                p = Graphics()
            Fxl = []
            Fyl = []
            Fxr = []
            Fyr = []    
            for i_a in range(n):
                for i_q in range(n):
                    sign_a = sign_l[i_a]()
                    sign_q = sign_r[i_q]()
                    if to_animate:
                        if sign_a > 0:
                            color_a = "red"
                        else:
                            color_a = "blue"
                        if sign_q > 0:
                            color_q = "red"
                        else:
                            color_q = "blue"

                    XYZa = p_l[i_a](t_i)
                    XYZq = p_r[i_q](t_i)
                    Xa = XYZa[0]
                    Ya = XYZa[1]
                    Za = XYZa[2]
                    Xq = XYZq[0]
                    Yq = XYZq[1]
                    Zq = XYZq[2]
                    if to_animate:
                        p += circle((cr[0],cr[1]), R_r)
                        p += circle((cl[0],cl[1]), R_l)
                        p += line ([(cr[0],cr[1]), (Xq, Yq)], color = color_q, linestyle="solid")
                        p += line ([(cl[0],cl[1]), (Xa, Ya)], color = color_a, linestyle="solid")
                    #print(Xa,Ya,Za)
                    (E_x, E_y, E_z, x_, y_, z_) = electr(Xa, Ya, Za, t_i, p_r[i_q], v_r[i_q], a_r[i_q], tlag=tlag_symbolic)
                    #print(E)
                    fxl = E_x*sign_a*sign_q
                    fyl = E_y*sign_a*sign_q
                    Fxl += [fxl]
                    Fyl += [fyl]
                    X2_q = x_
                    Y2_q = y_
                    (E_x, E_y, E_z, x_, y_, z_) = electr(Xq, Yq, Zq, t_i, p_l[i_a], v_l[i_a], a_l[i_a], tlag=tlag_symbolic)
                    #print(E)
                    fxr = E_x*sign_a*sign_q
                    fyr = E_y*sign_a*sign_q
                    Fxr += [fxr]
                    Fyr += [fyr]
                    X2_a = x_
                    Y2_a = y_
                    if to_animate:
                        p += line ([(Xq, Yq), (X2_a, Y2_a)], color = "green", linestyle="dashed")
                        p += line ([(X2_q, Y2_q), (Xa, Ya)], color = "green", linestyle="dashed")
                        p += line ([(cr[0],cr[1]), (X2_q, Y2_q)], color = color_q, linestyle="dashed")
                        p += line ([(cl[0],cl[1]), (X2_a, Y2_a)], color = color_a, linestyle="dashed")
            if to_animate:
                #p.show(aspect_ratio = 1)
                frames += [p]
            Fxt_l += [(t_i,Fxl)]
            Fyt_l += [(t_i,Fyl)]
            Fxt_r += [(t_i,Fxr)]
            Fyt_r += [(t_i,Fyr)]

    # Интегральная величина тяги в направлении оси y
    # угловое усилие

    return sum(Fx_r) + sum(Fx_l), sum(Fy_r) + sum(Fy_l), sum(F_alpha_l), sum(F_alpha_r), frames, Fxt_l, Fyt_l, Fxt_r, Fyt_r



# In[68]:

N = 16
for iA in range(0, N+1):
    Ai = iA * 2*pi / N
    print(Ai)


# In[ ]:


S_sum_Fx = []
S_sum_Fy = []
S_sum_F_alpha_l = []
S_sum_F_alpha_r = []
for iA in range(0, N):
    Ai = iA * 2*pi / N
    F_alpha = cals_sum_F(n=1, S=S, R_l=R_l, R_r=R_r, alpha0_l=Ai, to_animate = False)
    print(Ai, F_alpha[0], F_alpha[1], F_alpha[2], F_alpha[3])
    S_sum_Fx += [(Ai, F_alpha[0])]
    S_sum_Fy += [(Ai, F_alpha[1])]
    S_sum_F_alpha_l += [(Ai, F_alpha[2])]
    S_sum_F_alpha_r += [(Ai, F_alpha[3])]


# In[ ]:


list_plot(S_sum_Fx).save("sum_Fx.png")
list_plot(S_sum_Fy).save("sum_Fy.png")


# In[ ]:


list_plot(S_sum_F_alpha_l).save("F_alpha0_l.png")
list_plot(S_sum_F_alpha_r).save("F_alpha0_r.png")


# In[ ]:


F_x = num_int(lambda alpha0_i : cals_sum_F(n=1, S=S, R_l=R_l, R_r=R_r, alpha0_l=alpha0_i, to_animate = False)[0], 0, 2*pi)
print("F_x =", F_x)


# In[ ]:


F_y = num_int(lambda alpha0_i : cals_sum_F(n=1, S=S, R_l=R_l, R_r=R_r, alpha0_l=alpha0_i, to_animate = False)[1], 0, 2*pi)
print("F_y =", F_y)


# In[ ]:


# Definition for Polarization P and Magnetization M Fully Consistent with Maxwell’s Equations
# http://www.jpier.org/PIERB/pierb64/06.15100606.pdf

# https://www.stanfordmagnets.com/magnetization-of-rare-earth-magnets.html

# A good bar magnet may have a magnetic moment of magnitude 0.1 Am2 and a volume of 1 cm3, or 0.000001 m3, 
# and therefore an average magnetization magnitude is 100,000 A/m. 
# Iron can have a magnetization of around a million A/m. 
# Such a large value explains why magnets are so effective at producing magnetic fields.

# https://www.booksite.ru/fulltext/1/001/008/079/910.htm

# Намагниченность, характеристика магнитного состояния макроскопического физического тела; 
# в случае однородно намагниченного тела Н. определяется как магнитный момент J единицы объёма тела: 
# J = M/V, где М — магнитный момент тела, V — его объём. В случае неоднородно намагниченного тела Н. 
# определяется для каждой точки тела (точнее, для каждого физически малого объёма dV): J = dM/dV, 
# где dM — магнитный момент объёма dV. 

# Единица Н. в Международной системе единиц — ампер на метр (1 а/м — Н., 
# при которой 1 м^3 вещества обладает магнитным моментом 1 а×м^2), 
# в СГС системе единиц — эрг/(гс×см^3); 1 эрг/(гс×см^3) = 10^3 а/м.


# In[ ]:




