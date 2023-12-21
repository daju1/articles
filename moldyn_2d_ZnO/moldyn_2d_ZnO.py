from math import sqrt as sqrt
from math import sin as sin
from math import pi as pi

# ZnO	a = 3.25 c = 5.2	Wurtzite (HCP)
# ZnO	4.580	Halite (FCC)

l = 0.458 # nm
k = 10000
alpha = 5000
beta = 5000

class moldyn():

    def __init__(self, rows, cols):
        self.data = []
        self.rows = rows
        self.cols = cols
        self.atom_count = self.rows * self.cols

        self.tstep1 = 0.5;                      # [1.0e-15 s]
        self.tstep2 = self.tstep1 * self.tstep1;# [1.0e-30 s^2]

        self.vel = [];# [1.0e+3 m/s]
        self.acc = [];# [1.0e+12 m/s^2]

        self.mass = [];# [kg/mol]

        self.f = [];
        self.w = [];# workaround for bond's potential energy

        self.f_left_boundary = []
        self.acc_left_boundary = []
        self.w_left_boundary = []
        self.f_right_boundary = []
        self.acc_right_boundary = []
        self.w_right_boundary = []

        self.dkin_left_boundary = [0.0, 0.0]
        self.dkin_right_boundary = [0.0, 0.0]
        self.kin_left_boundary = [0.0, 0.0]
        self.kin_right_boundary = [0.0, 0.0]

        self.crd0 = [];
        self.crd = [];
        self.rc = [];
        self.nbr = [];
        self.nbri = [];

        self.step_counter = 0;

        self.sum_of_masses = 0.0;# [kg/mol]

        self.switch_xy = True

        counter = 0

        for row in range(self.rows):
            y0 = (row//4)*3*l
            if (row%4) == 0:
                y = 0 + y0
                x0 = 0
            elif (row%4) == 1:
                y = l/2 + y0
                x0 = l * (1/2*sqrt(3))
            elif (row%4) == 2:
                y = 3*l/2 + y0
                x0 = l * (1/2*sqrt(3))
            elif (row%4) == 3:
                y = 2*l + y0
                x0 = 0

            for col in range(self.cols):
                self.mass += [65.38 if row%2 == 0 else 15.999]
                self.mass[counter] *= 1.6605402e-27 * 6.0221367e+23;

                self.sum_of_masses += self.mass[counter];# kg/mol ; all atoms

                x = x0 + col * l * sqrt(3)

                self.rc.append([row, col])
                self.crd.append([x, -y])
                self.crd0.append([x, -y])

                #print(row, col, counter, x, y)

                self.vel.append([0.0, 0.0]);
                self.acc.append([0.0, 0.0]);

                self.f.append([0.0, 0.0]);
                self.w.append(0.0);

                self.f_left_boundary.append([0.0, 0.0]);
                self.acc_left_boundary.append([0.0, 0.0]);
                self.w_left_boundary.append(0.0);
                self.f_right_boundary.append([0.0, 0.0]);
                self.acc_right_boundary.append([0.0, 0.0]);
                self.w_right_boundary.append(0.0);

                if (row%4) == 0:
                    nb1 = (row - 1, col)
                    nb2 = (row + 1, col)
                    nb3 = (row + 1, col - 1)
                elif (row%4) == 1:
                    nb1 = (row - 1, col)
                    nb2 = (row - 1, col + 1)
                    nb3 = (row + 1, col)
                elif (row%4) == 2:
                    nb1 = (row - 1, col)
                    nb2 = (row + 1, col)
                    nb3 = (row + 1, col + 1)
                elif (row%4) == 3:
                    nb1 = (row - 1, col - 1)
                    nb2 = (row - 1, col)
                    nb3 = (row + 1, col)

                self.nbr.append([nb1, nb2, nb3])
                self.nbri.append([self.nbr2index(nb1), self.nbr2index(nb2), self.nbr2index(nb3)])

                counter+=1;

    def init_parabolic_dx(self, w, m):
        rw = self.rows // 2
        a = rw//2-(1+w)+rw%2
        b = rw//2+w

        for r in range(self.rows):
            i = r // 2
            dx = 0

            if i > a and i < b:
                dx = (i - a) * (b - i) / 16

            for c in range(self.cols):
                index = self.rc2index(r, c)
                self.crd[index][0] += l*dx*m

    def init_parabolic_dy(self, w, m):
        rw = self.rows // 2
        a = rw//2-(1+w)+rw%2
        b = rw//2+w

        for r in range(self.rows):
            i = r // 2
            dy = 0

            if i > a and i < b:
                dy = (i - a) * (b - i) / 16

            for c in range(self.cols):
                index = self.rc2index(r, c)
                self.crd[index][1] += l*dy * m

    def init_sin_dx(self, w, m):
        rw = self.rows // 2
        # w = 13
        # 2w = 26 -> 2*pi * 4 = 6.28 * 4
        a = rw//2 - w
        b = rw//2 + w
        T = w / pi
        center = rw//2
        print(a, b, center)

        for r in range(self.rows):
            i = r // 2
            dx = 0
            if i > a and i < b:
                dx = sin((i - center) / T)

            for c in range(self.cols):
                index = self.rc2index(r, c)
                self.crd[index][0] += l*dx*m

    def init_sin_dy(self, w, m):
        rw = self.rows // 2
        a = rw//2 - w
        b = rw//2 + w
        T = w / pi
        center = rw//2
        print(a, b, center)

        for r in range(self.rows):
            i = r // 2
            dy = 0
            if i > a and i < b:
                dy = sin((i - center) / T)

            for c in range(self.cols):
                index = self.rc2index(r, c)
                self.crd[index][1] += l*dy*m

    def init_sin_vx(self, w, m):
        rw = self.rows // 2
        # w = 13
        # 2w = 26 -> 2*pi * 4 = 6.28 * 4
        a = rw//2 - w
        b = rw//2 + w
        T = w / pi
        center = rw//2
        print(a, b, center)

        for r in range(self.rows):
            i = r // 2
            vx = 0
            if i > a and i < b:
                vx = sin((i - center) / T)

            for c in range(self.cols):
                index = self.rc2index(r, c)
                self.vel[index][0] += vx*m

    def init_sin_vy(self, w, m):
        rw = self.rows // 2
        # w = 13
        # 2w = 26 -> 2*pi * 4 = 6.28 * 4
        a = rw//2 - w
        b = rw//2 + w
        T = w / pi
        center = rw//2
        print(a, b, center)

        for r in range(self.rows):
            i = r // 2
            vy = 0
            if i > a and i < b:
                vy = sin((i - center) / T)

            for c in range(self.cols):
                index = self.rc2index(r, c)
                self.vel[index][1] += vy*m

    def nbr2index(self, nb):
        r = nb[0]
        c = nb[1]
        return self.rc2index(r, c)

    def rc2index(self, r, c):

        if r < 0: r = self.rows - 1
        if c < 0: c = self.cols - 1
        if r == self.rows: r = 0
        if c == self.cols: c = 0

        index = self.cols * r + c

        return index

    def get_zink_crd(self):
        if self.switch_xy:
            return [[self.crd[self.rc2index(r, c)][1], self.crd[self.rc2index(r, c)][0]] for r in range(self.rows) for c in range(self.cols) if r%2 == 0 ]
        return [self.crd[self.rc2index(r, c)] for r in range(self.rows) for c in range(self.cols) if r%2 == 0 ]

    def get_oxigen_crd(self):
        if self.switch_xy:
            return [[self.crd[self.rc2index(r, c)][1], self.crd[self.rc2index(r, c)][0]] for r in range(self.rows) for c in range(self.cols) if r%2 == 1 ]
        return [self.crd[self.rc2index(r, c)] for r in range(self.rows) for c in range(self.cols) if r%2 == 1 ]

    def box_x(self):
        return self.cols * l * sqrt(3)

    def box_y(self):
        return (self.rows//4) * 3 * l

    def boundary_x(self, crd1, crd2):
        if crd2[0] - crd1[0] > self.box_x() / 2:
            return -1
        if crd2[0] - crd1[0] < - self.box_x() / 2:
            return +1
        return 0

    def boundary_y(self, crd1, crd2):
        if crd2[1] - crd1[1] > self.box_y() / 2:
            return -1
        if crd2[1] - crd1[1] < - self.box_y() / 2:
            return +1
        return 0

    def distance_x(self, crd1, crd2):
        if crd2[0] - crd1[0] > self.box_x() / 2:
            return crd2[0] - crd1[0] - self.box_x()
        if crd2[0] - crd1[0] < - self.box_x() / 2:
            return crd2[0] - crd1[0] + self.box_x()
        return crd2[0] - crd1[0]

    def distance_y(self, crd1, crd2):
        if crd2[1] - crd1[1] > self.box_y() / 2:
            return crd2[1] - crd1[1] - self.box_y()
        if crd2[1] - crd1[1] < - self.box_y() / 2:
            return crd2[1] - crd1[1] + self.box_y()
        return crd2[1] - crd1[1]

    def distance(self, crd1, crd2):
        return sqrt(self.distance_x(crd1, crd2)**2 + self.distance_y(crd1, crd2)**2)

    def calc_atom(self, i):
        fx = 0.0
        fy = 0.0
        w = 0.0

        fy_left_boundary = 0.0
        wy_left_boundary = 0.0
        fy_right_boundary = 0.0
        wy_right_boundary = 0.0

        for nbr in self.nbri[i]:
            d  = self.distance  (self.crd[i], self.crd[nbr])
            dx = self.distance_x(self.crd[i], self.crd[nbr])
            dy = self.distance_y(self.crd[i], self.crd[nbr])
            
            bx = self.boundary_x(self.crd[i], self.crd[nbr])
            by = self.boundary_y(self.crd[i], self.crd[nbr])
            
            cx = dx/d
            cy = dy/d

            dl = d - l

            f = k * dl + alpha * dl**2 + beta * dl**3
            # bond's potential energy
            dw = k * dl * dl / 2 + alpha * dl**3 / 3 + beta * dl**4 / 4

            fx += f * cx
            fy += f * cy
            w  += dw

            if by > 0:
                fy_left_boundary += fy
                wy_left_boundary += dw/2
            if by < 0:
                fy_right_boundary += fy
                wy_right_boundary += dw/2

        self.f[i][0] = fx
        self.f[i][1] = fy
        self.w[i]    = w/2 # assign to atom half of bond's potential energy

        self.f_left_boundary[i][1]  = fy_left_boundary
        self.w_left_boundary[i]     = wy_left_boundary
        self.f_right_boundary[i][1] = fy_right_boundary
        self.w_right_boundary[i]    = wy_right_boundary

    def ComputeForce(self):
        for n1 in range(self.atom_count):
            self.calc_atom(n1)

    def CalcHeatFlow(self):
        for n2 in [0,1]:
            self.kin_left_boundary[n2]  += self.dkin_left_boundary[n2]
            self.kin_right_boundary[n2] += self.dkin_right_boundary[n2]
            self.dkin_left_boundary[n2]  = 0.0
            self.dkin_right_boundary[n2] = 0.0
            for row in [0, self.rows-1]:
                for col in range(self.cols):
                    n1 = self.rc2index(row, col)
                    tmpX = 500.0 * self.mass[n1]

                    acc_n1_n2 = self.f[n1][n2] / self.mass[n1];

                    self.acc_left_boundary[n1][n2]  = self.f_left_boundary[n1][n2] / self.mass[n1];
                    self.acc_right_boundary[n1][n2] = self.f_right_boundary[n1][n2] / self.mass[n1];
                    
                    # current velocity
                    cur_vel = self.vel[n1][n2]
                    # current partial kinetic energy
                    cur_kin = tmpX * cur_vel * cur_vel
                    
                    # updated velocity
                    new_vel = cur_vel + self.tstep1 * self.acc[n1][n2] * 1.0e-6
                    new_kin = tmpX * new_vel * new_vel
                    
                    # velocity updated just by left boundary force
                    new_vel_left_bound = cur_vel + self.tstep1 * self.acc_left_boundary[n1][n2] * 1.0e-6
                    new_kin_left_bound = tmpX * new_vel_left_bound * new_vel_left_bound
                    
                    # velocity updated just by right boundary force
                    new_vel_right_bound = cur_vel + self.tstep1 * self.acc_right_boundary[n1][n2] * 1.0e-6
                    new_kin_right_bound = tmpX * new_vel_right_bound * new_vel_right_bound
                    
                    self.dkin_left_boundary[n2]  += (new_kin_left_bound - new_kin) / self.tstep1 * 1.0e+15
                    self.dkin_right_boundary[n2] += (new_kin_right_bound - new_kin) / self.tstep1 * 1.0e+15
 
    def TakeMDStep(self):
        for n1 in range(self.atom_count):
            for n2 in [0,1]:
                tmpA = self.acc[n1][n2];

                tmp1 = self.tstep1 * self.vel[n1][n2] * 1.0e-3;
                tmp2 = self.tstep2 * tmpA * 0.5e-9;

                self.crd[n1][n2] += tmp1 + tmp2;
                self.vel[n1][n2] += self.tstep1 * tmpA * 0.5e-6;

        self.ComputeForce()
        self.CalcHeatFlow()

        for n1 in range(self.atom_count):
            for n2 in [0,1]:
                self.acc[n1][n2] = self.f[n1][n2] / self.mass[n1];
                self.vel[n1][n2] += self.tstep1 * self.acc[n1][n2] * 0.5e-6;

    def KineticEnergy(self):
        energy = 0.0

        for n1 in range(self.atom_count):
            tmpX = 500.0 * self.mass[n1]

            for n2 in [0,1]:
                tmp1 = self.vel[n1][n2];
                tmp2 = tmpX * tmp1 * tmp1;

                energy += tmp2;

        return energy;

    def KineticEnergy2(self):
        counter = 0
        energy1 = 0.0
        energy2 = 0.0
        for row in range(self.rows):
            for col in range(self.cols):
                tmpX = 500.0 * self.mass[counter]

                for n2 in [0,1]:
                    tmp1 = self.vel[counter][n2];
                    tmp2 = tmpX * tmp1 * tmp1;

                    if row < self.rows//2:
                        energy1 += tmp2;
                    else:
                        energy2 += tmp2;

                counter+=1;

        return energy1, energy2;

    def KineticEnergyCenter(self):
        counter = 0
        energy_center = [0.0, 0.0]

        for row in range(self.rows):
            for col in range(self.cols):
                tmpX = 500.0 * self.mass[counter]

                energy = 0.0
                for n2 in [0,1]:
                    tmp1 = self.vel[counter][n2];
                    tmp2 = tmpX * tmp1 * tmp1;
                    energy += tmp2
                    

                for n2 in [0,1]:
                    crd = self.crd[counter][n2]
                    energy_center[n2] += crd * energy

                counter+=1;

        for n2 in [0,1]:
            energy_center[n2] /= self.KineticEnergy()

        return energy_center

    def PotentialEnergy(self):
        counter = 0
        energy = 0.0
        for row in range(self.rows):
            for col in range(self.cols):
                tmpW = self.w[counter]
                energy += tmpW;
                counter+=1;

        return energy

    def PotentialEnergy2(self):
        counter = 0
        energy1 = 0.0
        energy2 = 0.0
        for row in range(self.rows):
            for col in range(self.cols):
                tmpW = self.w[counter]

                if row < self.rows//2:
                    energy1 += tmpW;
                else:
                    energy2 += tmpW;

                counter+=1;

        return energy1, energy2

    def PotentialEnergyCenter(self):
        counter = 0
        energy_center = [0.0, 0.0]

        for row in range(self.rows):
            for col in range(self.cols):
                tmpW = self.w[counter]

                for n2 in [0,1]:
                    crd = self.crd[counter][n2]

                    energy_center[n2] += crd * tmpW

                counter+=1;

        for n2 in [0,1]:
            energy_center[n2] /= self.PotentialEnergy()

        return energy_center

    def FullEnergyCenter(self):
        counter = 0
        energy_center = [0.0, 0.0]

        for row in range(self.rows):
            for col in range(self.cols):
                tmpX = 500.0 * self.mass[counter]
                tmpW = self.w[counter]

                energy = 0.0
                for n2 in [0,1]:
                    tmp1 = self.vel[counter][n2];
                    tmp2 = tmpX * tmp1 * tmp1;
                    energy += tmp2

                for n2 in [0,1]:
                    crd = self.crd[counter][n2]
                    energy_center[n2] += crd * (energy + tmpW)

                counter+=1;

        for n2 in [0,1]:
            energy_center[n2] /= (self.KineticEnergy()+self.PotentialEnergy())

        return energy_center


    def LagrangianDencityCenter(self):
        counter = 0
        energy_center = [0.0, 0.0]

        for row in range(self.rows):
            for col in range(self.cols):
                tmpX = 500.0 * self.mass[counter]
                tmpW = self.w[counter]

                energy = 0.0
                for n2 in [0,1]:
                    tmp1 = self.vel[counter][n2];
                    tmp2 = tmpX * tmp1 * tmp1;
                    energy += tmp2

                for n2 in [0,1]:
                    crd = self.crd[counter][n2]
                    # L = T - U
                    energy_center[n2] += crd * (energy - tmpW)

                counter+=1;
                
        for n2 in [0,1]:
            energy_center[n2] /= (self.KineticEnergy()-self.PotentialEnergy())
            
        return energy_center

    def plot(self):
        from sage.plot.plot import Graphics
        from sage.plot.point import point
        p = Graphics()
        p += point(self.get_zink_crd(),   marker='o', markeredgecolor='blue', size=20)
        p += point(self.get_oxigen_crd(), marker='o', markeredgecolor='red',  size=2)
        return p

    def plot_dx(self, dx_min, dx_max):
        from sage.plot.plot import list_plot
        plot_data = []
        c = self.cols // 2
        for r in range(self.rows):
            index = self.rc2index(r, c)
            x0 = self.crd0[index][0]
            y0 = self.crd0[index][1]
            dx = self.crd[index][0] - x0
            plot_data += [(-y0, dx)]

            if dx_min > dx:
                dx_min = dx
            if dx_max < dx:
                dx_max = dx

        p = list_plot(plot_data, size=20, axes_labels=['$y$','$dx$'])
        return p, dx_min, dx_max

    def plot_dy(self, dy_min, dy_max):
        from sage.plot.plot import list_plot
        plot_data = []
        c = self.cols // 2
        for r in range(self.rows):
            index = self.rc2index(r, c)
            y0 = self.crd0[index][1]
            dy = self.crd[index][1] - y0
            plot_data += [(-y0, dy)]

            if dy_min > dy:
                dy_min = dy
            if dy_max < dy:
                dy_max = dy

        p = list_plot(plot_data, size=20, axes_labels=['$y$','$dy$'])
        return p, dy_min, dy_max

    def plot_vx(self, vx_min, vx_max):
        from sage.plot.plot import list_plot
        plot_data = []
        c = self.cols // 2
        for r in range(self.rows):
            index = self.rc2index(r, c)
            y0 = self.crd0[index][1]
            vx = self.vel[index][0]
            plot_data += [(-y0, vx)]

            if vx_min > vx:
                vx_min = vx
            if vx_max < vx:
                vx_max = vx

        p = list_plot(plot_data, size=20, axes_labels=['$y$','$v_x$'])
        return p, vx_min, vx_max

    def plot_vy(self, vy_min, vy_max):
        from sage.plot.plot import list_plot
        plot_data = []
        c = self.cols // 2
        for r in range(self.rows):
            index = self.rc2index(r, c)
            y0 = self.crd0[index][1]
            vy = self.vel[index][1]
            plot_data += [(-y0, vy)]

            if vy_min > vy:
                vy_min = vy
            if vy_max < vy:
                vy_max = vy

        p = list_plot(plot_data, size=20, axes_labels=['$y$','$v_y$'])
        return p, vy_min, vy_max

    def xmin(self):
        return -l/2 * (1/2*sqrt(3))
    def xmax(self):
        return -l/2 * (1/2*sqrt(3)) + self.box_x()

    def ymax(self):
        return +l/2
    def ymin(self):
        return l/2 - self.box_y()

    def show(self, p):
        if self.switch_xy:
            p.show(aspect_ratio=1, axes=False,
                             ymax = self.xmax(), ymin = self.xmin(),
                             xmin = self.ymin(), xmax = self.ymax()
                            )
        else:
            p.show(aspect_ratio=1, axes=False,
                             ymax = self.ymax(), ymin = self.ymin(),
                             xmin = self.xmin(), xmax = self.xmax()
                            )