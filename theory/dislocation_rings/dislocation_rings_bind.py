import ctypes, math
lib = ctypes.CDLL("dislocation_rings.so")

class RingConfig(ctypes.Structure):
    _fields_ = [("center", ctypes.c_double*3),
                ("e1", ctypes.c_double*3),
                ("e2", ctypes.c_double*3),
                ("radius", ctypes.c_double),
                ("omega", ctypes.c_double),
                ("phase0", ctypes.c_double),
                ("total_charge", ctypes.c_double),
                ("n_elem", ctypes.c_int),
                ("charge_sign", ctypes.c_int)]

lib.compute_ring_dipole_force.argtypes = [
    ctypes.POINTER(RingConfig),
    ctypes.c_int,
    ctypes.POINTER(RingConfig),
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_int,
    ctypes.c_double*3
]

lib.compute_ring_dipole_force.restype = None

def make_ring(
        center,
        radius,
        omega,
        phase0,
        total_charge,
        n_elem=24,
        sign=1):
    r = RingConfig();
    r.center[:]=center;
    r.e1[:]=(1.0,0.0,0.0);
    r.e2[:]=(0.0,1.0,0.0)
    r.radius=radius;
    r.omega=omega;
    r.phase0=phase0;
    r.total_charge=total_charge;
    r.n_elem=n_elem;
    r.charge_sign=sign
    return r

def force_multi(
        field_rings,
        source_rings,
        t, c, mode):
    F=(ctypes.c_double*3)()
    FR=(RingConfig*len(field_rings))(*field_rings)
    SR=(RingConfig*len(source_rings))(*source_rings)
    lib.compute_ring_dipole_force(
        FR,len(field_rings),
        SR,len(source_rings),
        t,c,mode,F)
    return [F[0],F[1],F[2]]

c=1.0;
a=0.5;
omega=0.3;
d=0.3

def make_dipole(
        center,
        dphi=0.0,
        sign_north=1):

    ring_N = make_ring(
        (center[0],center[1],center[2]+d/2),
        a,
        omega,
        dphi,
        1.0,
        sign=sign_north
    )
    ring_S = make_ring(
        (center[0],center[1],center[2]-d/2),
        a,
        omega,
        dphi,
        1.0,
        sign=-sign_north)

    return [ring_N, ring_S]
