from sympy import Point, Line, Segment, Ray
import numpy as np

import sage
from sage.symbolic.constants import *
from sage.plot.graphics import *
from sage.plot.line import line
from sage.plot.arrow import arrow

from enum import Enum

class State(Enum):
    STARTED = 2
    INTERSECTED_DOWN = 3
    INTERSECTED_UPPER = 4
    CAME_OUT = 1
    CAME_IN = 0


epsilon = np.float128(1e-16)

color = "green"

calc_float = True
draw_normal = False
logging = False


def draw_wedge(alpha, r, d, traj_ray_1, draw_plot=True):
    # d толщина мембраны
    # r половина размера меньшей щели, которая слева
    # alpha образующий угол клина
    a = np.float128(0.1) # дополнительный размер плоскости мембраны со стороны большего размера отверстия

    if draw_plot:
        plt = sage.plot.graphics.Graphics()
    else:
        plt = None

    # R половина размера большей щели, которая справа
    R = np.float128((r + d * np.tan(np.float128(alpha))))

    p1 = (0, r)
    p2 = (d, R)
    wedge_segment_upper = Segment(p1, p2)
    if draw_plot:
        plt += line (wedge_segment_upper.points, color = color, linestyle="dashed")

    p1 = (0, r)
    p2 = (0, R + a)
    if draw_plot:
        plt += line ([p1, p2], color = color, linestyle="dashed")

    p1 = (d, R)
    p2 = (d, R + a)
    if draw_plot:
        plt += line ([p1, p2], color = color, linestyle="dashed")



    p1 = (0, -r)
    p2 = (d, -R)
    wedge_segment_down = Segment(p1, p2)
    if draw_plot:
        plt += line (wedge_segment_down.points, color = color, linestyle="dashed")

    p1 = (0, -r)
    p2 = (0, -R - a)
    if draw_plot:
        plt += line ([p1, p2], color = color, linestyle="dashed")

    p1 = (d, -R)
    p2 = (d, -R - a)
    if draw_plot:
        plt += line ([p1, p2], color = color, linestyle="dashed")


    segment_out = Segment((np.float128(0),  r),
                          (np.float128(0), -r))
    segment_in  = Segment((d,  R),
                          (d, -R))

    state = State.STARTED
    number_of_intersections = 0

    def traj_ray_intersect_wedge_segment_down(plt, alpha, traj_ray, number_of_intersections):
        if logging:
            print("traj_ray_intersect_wedge_segment_down")
        try:
            intersection_down = traj_ray.intersect(wedge_segment_down)
            if logging:
                print("intersection_down",
                      intersection_down)
        except Exception as ex:
            print (ex)
            intersection_down = ()

        if (intersection_down == wedge_segment_down):
            if traj_ray.xdirection > 0:
                return plt, State.CAME_IN, number_of_intersections
            elif traj_ray.xdirection < 0:
                return plt, State.CAME_OUT, number_of_intersections
        elif len(intersection_down) > 0:
            intersection_point_down, = intersection_down

            if logging:
                print("we got intersection_point_down",
                      intersection_point_down)

            #if traj_ray.source != intersection_point_down:
            if  abs(traj_ray.source.x - intersection_point_down.x) > epsilon or \
                abs(traj_ray.source.y - intersection_point_down.y) > epsilon:

                if calc_float:
                    intersection_point_down = Point(
                        np.float128(intersection_point_down.x),
                        np.float128(intersection_point_down.y))

                if logging:
                    print("intersection_point_down",
                          intersection_point_down)

                if draw_plot:
                    plt += arrow(traj_ray.source, intersection_point_down, linestyle="dashed")
                    #print("arrow", traj_ray.source, intersection_point_down)

                traj_ray_angle_between_wedge_segment_down = traj_ray.angle_between(wedge_segment_down)
                if calc_float:
                    traj_ray_angle_between_wedge_segment_down = np.float128(
                        traj_ray_angle_between_wedge_segment_down)

                if logging:
                    print("traj_ray.angle_between(wedge_segment_down)",
                          traj_ray_angle_between_wedge_segment_down)

                down_normal_angle = pi/2-alpha
                if logging:
                    print("down_normal_angle", (down_normal_angle/pi).n())

                if draw_normal:
                    ray_normal = Ray(intersection_point_down, angle = down_normal_angle)
                    if logging:
                        print("ray_normal"
                              , np.float128(ray_normal.points[0].x)
                              , np.float128(ray_normal.points[0].y)
                              , np.float128(ray_normal.points[1].x)
                              , np.float128(ray_normal.points[1].y)
                            )

                    intersection_normal = ray_normal.intersect(wedge_segment_upper)
                    if len(intersection_normal) > 0:
                        intersection_normal_point, = intersection_normal

                        intersection_normal_point = Point(
                            np.float128(intersection_normal_point.x),
                            np.float128(intersection_normal_point.y))
                        if logging:
                            print("intersection_normal_point with wedge_segment_upper", intersection_normal_point)
                        if draw_plot:
                            plt += line((ray_normal.source, intersection_normal_point))
                    else:
                        intersection_normal = ray_normal.intersect(segment_in)
                        if len(intersection_normal) > 0:
                            intersection_normal_point, = intersection_normal

                            intersection_normal_point = Point(
                                np.float128(intersection_normal_point.x),
                                np.float128(intersection_normal_point.y))
                            if logging:
                                print("intersection_normal_point with wedge_segment_upper", intersection_normal_point)
                            if draw_plot:
                                plt += line((ray_normal.source, intersection_normal_point))

                angle_ray_2 = down_normal_angle + (traj_ray_angle_between_wedge_segment_down - pi/2)
                if calc_float:
                    angle_ray_2 = np.float128(angle_ray_2)
                if logging:
                    print("angle_ray_2 from down", (angle_ray_2/pi).n())
                traj_ray_2 = Ray(intersection_point_down, angle = angle_ray_2)

                number_of_intersections += 1

                return draw_traj(plt, alpha, traj_ray_2, State.INTERSECTED_DOWN, number_of_intersections)

            #if traj_ray.source == intersection_point_down:
            else:
                traj_ray_angle_between_wedge_segment_down = traj_ray.angle_between(wedge_segment_down)
                if calc_float:
                    traj_ray_angle_between_wedge_segment_down = np.float128(
                        traj_ray_angle_between_wedge_segment_down)

                if logging:
                    print("traj_ray.angle_between(wedge_segment_down)",
                          traj_ray_angle_between_wedge_segment_down)

                down_normal_angle = pi/2-alpha
                if logging:
                    print("down_normal_angle", down_normal_angle)
                angle_ray_2 = down_normal_angle + \
                    (traj_ray_angle_between_wedge_segment_down - pi/2)
                if calc_float:
                    angle_ray_2 = np.float128(angle_ray_2)
                if logging:
                    print("angle_ray_2", (angle_ray_2/pi).n())
                traj_ray_2 = Ray(traj_ray.source, angle = angle_ray_2)

                if logging:
                    if draw_plot:
                        plt += arrow(traj_ray_2.source,
                                 traj_ray_2.points[1],
                                 color="green",
                                 thickness = 0.1
                                )

                number_of_intersections += 1

                return draw_traj(plt, alpha, traj_ray_2, State.INTERSECTED_DOWN, number_of_intersections)

    def traj_ray_intersect_wedge_segment_upper(plt, alpha, traj_ray, number_of_intersections):
        if logging:
            print("traj_ray_intersect_wedge_segment_upper")
        try:
            intersection_upper = traj_ray.intersect(wedge_segment_upper)
            if logging:
                print("intersection_upper",
                      intersection_upper)
        except Exception as ex:
            print (ex)
            intersection_upper = ()

        if (intersection_upper == wedge_segment_upper):
            if traj_ray.xdirection > 0:
                return plt, State.CAME_IN, number_of_intersections
            elif traj_ray.xdirection < 0:
                return plt, State.CAME_OUT, number_of_intersections
        elif len(intersection_upper) > 0:
            intersection_point_upper, = intersection_upper

            if logging:
                print("we got intersection_point_upper",
                      intersection_point_upper)

            # if traj_ray.source != intersection_point_upper:
            if  abs(traj_ray.source.x - intersection_point_upper.x) > epsilon or \
                abs(traj_ray.source.y - intersection_point_upper.y) > epsilon:

                if calc_float:
                    intersection_point_upper = Point(
                        np.float128(intersection_point_upper.x),
                        np.float128(intersection_point_upper.y))
                if logging:
                    print("intersection_point_upper",
                          intersection_point_upper)

                if draw_plot:
                    plt += arrow(traj_ray.source, intersection_point_upper, linestyle="dashed")
                    #print("arrow", traj_ray.source, intersection_point_upper)


                traj_ray_angle_between_wedge_segment_upper = traj_ray.angle_between(wedge_segment_upper)
                if calc_float:
                    traj_ray_angle_between_wedge_segment_upper = np.float128(
                        traj_ray_angle_between_wedge_segment_upper)
                if logging:
                    print("traj_ray.angle_between(wedge_segment_upper)",
                          traj_ray_angle_between_wedge_segment_upper)

                upper_normal_angle = 2*pi-pi/2+alpha
                if logging:
                    print("upper_normal_angle", (upper_normal_angle/pi).n())

                if draw_normal:
                    ray_normal = Ray(intersection_point_upper, angle = upper_normal_angle)
                    if logging:
                        print("ray_normal"
                          , np.float128(ray_normal.points[0].x)
                          , np.float128(ray_normal.points[0].y)
                          , np.float128(ray_normal.points[1].x)
                          , np.float128(ray_normal.points[1].y)
                        )

                    intersection_normal = ray_normal.intersect(wedge_segment_down)
                    if len(intersection_normal) > 0:
                        intersection_normal_point, = intersection_normal

                        intersection_normal_point = Point(
                            np.float128(intersection_normal_point.x),
                            np.float128(intersection_normal_point.y))
                        if logging:
                            print("intersection_normal_point with wedge_segment_down",
                                  intersection_normal_point)
                        if draw_plot:
                            plt += line((ray_normal.source, intersection_normal_point))
                    else:
                        intersection_normal = ray_normal.intersect(segment_in)
                        if len(intersection_normal) > 0:
                            intersection_normal_point, = intersection_normal

                            intersection_normal_point = Point(
                                np.float128(intersection_normal_point.x),
                                np.float128(intersection_normal_point.y))
                            if logging:
                                print("intersection_normal_point with wedge_segment_upper",
                                      intersection_normal_point)
                            if draw_plot:
                                plt += line((ray_normal.source, intersection_normal_point))

                angle_ray_2 = upper_normal_angle - (traj_ray_angle_between_wedge_segment_upper - pi/2)
                if calc_float:
                    angle_ray_2 = np.float128(angle_ray_2)
                if logging:
                    print("angle_ray_2 from upper", (angle_ray_2/pi).n())
                traj_ray_2 = Ray(intersection_point_upper, angle = angle_ray_2)

                number_of_intersections += 1

                return draw_traj(plt, alpha, traj_ray_2, State.INTERSECTED_UPPER, number_of_intersections)

            # if traj_ray.source == intersection_point_upper:
            else:
                traj_ray_angle_between_wedge_segment_upper = traj_ray.angle_between(wedge_segment_upper)
                if calc_float:
                    traj_ray_angle_between_wedge_segment_upper = np.float128(
                        traj_ray_angle_between_wedge_segment_upper)
                if logging:
                    print("traj_ray.angle_between(wedge_segment_upper)",
                          traj_ray_angle_between_wedge_segment_upper)

                upper_normal_angle = 2*pi-pi/2+alpha
                if logging:
                    print("upper_normal_angle", upper_normal_angle)
                angle_ray_2 = upper_normal_angle - \
                    (traj_ray_angle_between_wedge_segment_upper - pi/2)
                if calc_float:
                    angle_ray_2 = np.float128(angle_ray_2)
                if logging:
                    print("angle_ray_2", angle_ray_2)
                traj_ray_2 = Ray(traj_ray.source, angle = angle_ray_2)

                if logging:
                    if draw_plot:
                        plt += arrow(traj_ray_2.source,
                                 traj_ray_2.points[1],
                                 color="green",
                                 thickness = 0.1)

                number_of_intersections += 1

                return draw_traj(plt, alpha, traj_ray_2, State.INTERSECTED_UPPER, number_of_intersections)

    def traj_ray_intersect_segment_in(plt, alpha, traj_ray, number_of_intersections):
        try:
            intersection_in = traj_ray.intersect(segment_in)
        except:
            intersection_in = ()

        if intersection_in == segment_in:
            return plt, State.CAME_IN, number_of_intersections
        elif len(intersection_in) > 0:
            intersection_point_in, = intersection_in

            #if traj_ray.source != intersection_point_in:
            if  abs(traj_ray.source.x - intersection_point_in.x) > epsilon or \
                abs(traj_ray.source.y - intersection_point_in.y) > epsilon:

                if draw_plot:
                    plt += arrow(traj_ray.source, intersection_point_in, linestyle="dashed")
                    #print("arrow", traj_ray.source, intersection_point_in)


                return plt, State.CAME_IN, number_of_intersections
            else:
                if traj_ray.xdirection > 0:
                    return plt, State.CAME_IN, number_of_intersections
                if logging:
                    print("traj_ray.source", traj_ray.source)
                    print("intersection_point_in", intersection_point_in)

    def traj_ray_intersect_segment_out(plt, alpha, traj_ray, number_of_intersections):
        try:
            intersection_out = traj_ray.intersect(segment_out)
        except:
            intersection_out = ()

        if intersection_out == segment_out:
            return plt, State.CAME_OUT, number_of_intersections
        elif len(intersection_out) > 0:
            intersection_point_out, = intersection_out

            #if traj_ray.source != intersection_point_out:
            if  abs(traj_ray.source.x - intersection_point_out.x) > epsilon or \
                abs(traj_ray.source.y - intersection_point_out.y) > epsilon:

                if draw_plot:
                    plt += arrow(traj_ray.source, intersection_point_out, linestyle="dashed")
                    #print("arrow", traj_ray.source, intersection_point_out)

                return plt, State.CAME_OUT, number_of_intersections
            #if traj_ray.source == intersection_point_out:
            else:
                if traj_ray.xdirection < 0:
                    return plt, State.CAME_OUT, number_of_intersections
                if logging:
                    print("traj_ray.source", traj_ray.source)
                    print("intersection_point_out", intersection_point_out)

    def draw_traj(plt, alpha, traj_ray, state, number_of_intersections):
        if logging:
            if draw_plot:
                plt.show(aspect_ratio = 1, axes=False)
            print("traj_ray", traj_ray)
            print("traj_ray.slope", traj_ray.slope)
            print("traj_ray.direction", traj_ray.direction)
            print("traj_ray.xdirection", traj_ray.xdirection)
            print("traj_ray.ydirection", traj_ray.ydirection)
            print("wedge_segment_down.slope", wedge_segment_down.slope)
            print("wedge_segment_upper.slope", wedge_segment_upper.slope)
            print("segment_in.slope", segment_in.slope)
            print("segment_out.slope", segment_out.slope)

        if traj_ray.xdirection > 0:
            ret = traj_ray_intersect_segment_in(plt, alpha, traj_ray, number_of_intersections)
            if None != ret:
                return ret

        if traj_ray.xdirection < 0:
            ret = traj_ray_intersect_segment_out(plt, alpha, traj_ray, number_of_intersections)
            if None != ret:
                return ret

        if traj_ray.ydirection < 0 and traj_ray.xdirection <= 0:
            ret = traj_ray_intersect_wedge_segment_down(plt, alpha, traj_ray, number_of_intersections)
            if None != ret:
                return ret

        elif traj_ray.ydirection < 0 and traj_ray.xdirection > 0 \
        and wedge_segment_down.slope > traj_ray.slope:
            ret = traj_ray_intersect_wedge_segment_down(plt, alpha, traj_ray, number_of_intersections)
            if None != ret:
                return ret

        elif traj_ray.xdirection < 0 \
        and wedge_segment_down.slope < traj_ray.slope:
            ret = traj_ray_intersect_wedge_segment_down(plt, alpha, traj_ray, number_of_intersections)
            if None != ret:
                return ret


        if traj_ray.ydirection > 0 and traj_ray.xdirection <= 0:
            ret = traj_ray_intersect_wedge_segment_upper(plt, alpha, traj_ray, number_of_intersections)
            if None != ret:
                return ret

        elif traj_ray.ydirection > 0 and traj_ray.xdirection > 0 \
        and wedge_segment_upper.slope < traj_ray.slope:
            ret = traj_ray_intersect_wedge_segment_upper(plt, alpha, traj_ray, number_of_intersections)
            if None != ret:
                return ret

        elif traj_ray.xdirection < 0 \
        and wedge_segment_upper.slope > traj_ray.slope:
            ret = traj_ray_intersect_wedge_segment_upper(plt, alpha, traj_ray, number_of_intersections)
            if None != ret:
                return ret

        #if logging:
        if draw_plot:
            plt += arrow(traj_ray.source, traj_ray.points[1], color="red")
            plt.show(aspect_ratio = 1, axes=False)

        return None


    return draw_traj(plt, alpha, traj_ray_1, state, number_of_intersections)



def integrand_in (alpha, r, d, y, ang, draw_plot=True):
    calc_float = True
    return draw_wedge(np.float128(alpha), np.float128(r), np.float128(d), traj_ray_1 = Ray((np.float128(d), np.float128(y)), angle = np.float128(ang)), draw_plot=draw_plot)



def integrand_out (alpha, r, d, y, ang, draw_plot=True):
    calc_float = True
    return draw_wedge(np.float128(alpha), np.float128(r), np.float128(d), traj_ray_1 = Ray((np.float128(0), np.float128(y)), angle = np.float128(ang)), draw_plot=draw_plot)



def make_investigation_plots(alpha, r, d, milk_angle, draw_plot=True):
    import numpy as np

    R = np.float128((r + d * np.tan(np.float128(alpha))))

    y_linspace   = np.linspace(-R, R, 101)
    ang_linspace = np.linspace(np.float128(pi/2+milk_angle),
                               np.float128(pi+pi/2-milk_angle),
                               361)

    y_list = y_linspace.tolist()
    ang_list = ang_linspace.tolist()

    y_grid, ang_grid = np.meshgrid(y_linspace, ang_linspace)

    s = y_grid * np.nan
    n = y_grid * np.nan

    frames = []

    for iy in np.arange(0, len(y_linspace), 1):
        for iang in np.arange(0, len(ang_linspace), 1):
            y = y_list[iy]
            ang = ang_list[iang]
            plt, state, number_of_reflections = \
                integrand_in(alpha=alpha, r=r, d=d, y=y, ang=ang, draw_plot=draw_plot)

            if 0 == iy % 2 and 0 == iang % 4:
                frames += [plt]
            if State.CAME_OUT == state:
                s[iang][iy] = 1
            elif State.CAME_IN == state:
                s[iang][iy] = 0
            n[iang][iy] = number_of_reflections
    return s, n, y_grid, ang_grid, frames



def calc_passing_coefficient(s, alpha, r, d, milk_angle):

    R = np.float128((r + d * np.tan(np.float128(alpha))))

    sum_passed = 0
    sum_all = 0
    for iy in np.arange(0, len(s[0]), 1):
        for iang in np.arange(0, len(s), 1):
            if not np.isnan(s[iang][iy]):
                sum_passed += s[iang][iy]
                sum_all += 1

    # increase sum all because
    # I played angles not from pi/2 till pi+pi/2
    # but from pi/2+milk_angle till pi+pi/2-milk_angle
    sum_all *= pi / (pi - 2*milk_angle)

    print ("sum passed", sum_passed)
    print ("sum all", sum_all)
    print ("sum passed/all", sum_passed/sum_all)

    # how differs squares of wedge sided in / out
    s_in_per_out = R/r

    print ("relation of squares of wedge sided in / out", s_in_per_out)

    passing_coefficient = sum_passed/sum_all * s_in_per_out
    print ("passing coefficient", passing_coefficient)

    return passing_coefficient



def plot_state_depending_on_angle(y, alpha, r, d, milk_angle):
    state_y0 = []
    reflections_y0 = []
    R = np.float128((r + d * np.tan(np.float128(alpha))))

    ang_linspace = np.linspace(np.float128(pi/2+milk_angle),
                               np.float128(pi+pi/2-milk_angle),
                               361)

    ang_list = ang_linspace.tolist()

    state_v = 0

    for iang in np.arange(0, len(ang_linspace), 1):

        ang = ang_list[iang]
        try:
            plt, state, number_of_reflections = \
                integrand_in (alpha=alpha, r=r, d=d, y=y, ang=ang, draw_plot=False)

            if State.CAME_OUT == state:
                state_v = 1.1
            elif State.CAME_IN == state:
                state_v = 0

            state_y0 += [(ang, state_v)]
            reflections_y0 += [(ang, number_of_reflections)]

        except:
            pass
    return state_y0, reflections_y0




def iterate_search_switching_angle_dw(y, alpha, r, d, draw_plot=True):

    epsilon = pi / 4096 / 4096 / 4096

    angle_switch_dw = np.arctan2(-np.float128(y), -np.float128(d))
    angle_switch_dw_plus_delta = np.arctan2(-np.float128(y+r), -np.float128(d))
    delta_angle_dw = angle_switch_dw_plus_delta - angle_switch_dw
    if False:
        print("angle_switch_dw", angle_switch_dw)
        print("angle_switch_dw_plus_delta", angle_switch_dw_plus_delta)
        print("delta_angle_dw", delta_angle_dw)
    if delta_angle_dw < 0:
        delta_angle_dw += np.float128(2*pi)
        #print("delta_angle_dw", delta_angle_dw)
    state_dw = State.CAME_OUT

    while(True):

        ang = angle_switch_dw

        plt, state, number_of_reflections = \
            integrand_in (alpha=alpha, r=r, d=d, y=y, ang=ang, draw_plot=draw_plot)

        if State.CAME_OUT == state:
            if State.CAME_OUT == state_dw:
                pass
            elif State.CAME_IN == state_dw:
                delta_angle_dw /= 2
            angle_switch_dw += delta_angle_dw
            state_dw = state

        elif State.CAME_IN == state:
            if State.CAME_OUT == state_dw:
                delta_angle_dw /= 2
            elif State.CAME_IN == state_dw:
                pass
            angle_switch_dw -= delta_angle_dw
            state_dw = state

        if draw_plot:
            plt.show(aspect_ratio = 1, axes=False)
        if epsilon > delta_angle_dw:
            break

    return angle_switch_dw, state_dw, delta_angle_dw



def iterate_search_switching_angle_up(y, alpha, r, d, draw_plot=True):

    epsilon = pi / 4096 / 4096 / 4096

    angle_switch_up = np.arctan2(-np.float128(y), -np.float128(d))
    angle_switch_up_plus_delta = np.arctan2(-np.float128(y-r), -np.float128(d))
    delta_angle_up = angle_switch_up - angle_switch_up_plus_delta
    if False:
        print("angle_switch_up", angle_switch_up)
        print("angle_switch_up_plus_delta", angle_switch_up_plus_delta)
        print("delta_angle_up", delta_angle_up)

    if delta_angle_up < 0:
        delta_angle_up += np.float128(2*pi)
        #print("delta_angle_up", delta_angle_up)
    state_up = State.CAME_OUT

    while(True):

        ang = angle_switch_up

        plt, state, number_of_reflections = \
            integrand_in (alpha=alpha, r=r, d=d, y=y, ang=ang, draw_plot=draw_plot)

        if State.CAME_OUT == state:
            if State.CAME_OUT == state_up:
                pass
            elif State.CAME_IN == state_up:
                delta_angle_up /= 2
            angle_switch_up -= delta_angle_up
            state_up = state

        elif State.CAME_IN == state:
            if State.CAME_OUT == state_up:
                delta_angle_up /= 2
            elif State.CAME_IN == state_up:
                pass
            angle_switch_up += delta_angle_up
            state_up = state

        if draw_plot:
            plt.show(aspect_ratio = 1, axes=False)
        if epsilon > delta_angle_up:
            break

    return angle_switch_up, state_up, delta_angle_up


def get_iterate_passing_angle_depending_on_y(y, alpha, r, d, draw_plot=False):
    angle_switch_up, state_up, delta_angle_up = iterate_search_switching_angle_up(y=y, alpha=alpha, r=r, d=d, draw_plot=draw_plot)
    angle_switch_dw, state_dw, delta_angle_dw = iterate_search_switching_angle_dw(y=y, alpha=alpha, r=r, d=d, draw_plot=draw_plot)
    return angle_switch_dw - angle_switch_up
