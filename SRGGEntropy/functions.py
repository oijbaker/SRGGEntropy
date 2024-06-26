import numpy as np
import matplotlib.pyplot as plt
import ctypes
import os.path

dll_name = "clibrary.so"
dllabspath = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + dll_name
libObject = ctypes.CDLL(dllabspath)

def square_entropy_hard(n, r0, lim):
    sq_ent_hard = libObject.square_entropy_hard
    sq_ent_hard.restype = ctypes.c_double
    sq_ent_hard.argtypes = [ctypes.c_int, ctypes.c_double,ctypes.c_int]

    return sq_ent_hard(n, r0, lim)

def square_entropy_triangular(n, r0, lim):
    sq_ent_tri = libObject.square_entropy_triangular
    sq_ent_tri.restype = ctypes.c_double
    sq_ent_tri.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return sq_ent_tri(n, r0, lim)

def square_entropy_rayleigh(n, r0, eta, lim):
    sq_ent_ray = libObject.square_entropy_rayleigh
    sq_ent_ray.restype = ctypes.c_double
    sq_ent_ray.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]

    return sq_ent_ray(n, r0, eta, lim)

def square_entropy_exclusion(n, r0, lim):
    sq_ent_excl = libObject.square_entropy_exclusion
    sq_ent_excl.restype = ctypes.c_double
    sq_ent_excl.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return sq_ent_excl(n, r0, lim)


def torus_1d_entropy_rayleigh(n, r0, eta, lim):
    torus_1d_ent_ray = libObject.torus_1d_entropy_rayleigh
    torus_1d_ent_ray.restype = ctypes.c_double
    torus_1d_ent_ray.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]

    return torus_1d_ent_ray(n, r0, eta, lim)

def torus_1d_entropy_hard(n, r0, lim):
    torus_1d_ent_hard = libObject.torus_1d_entropy_hard
    torus_1d_ent_hard.restype = ctypes.c_double
    torus_1d_ent_hard.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return torus_1d_ent_hard(n, r0, lim)

def torus_1d_entropy_triangular(n, r0, lim):
    torus_1d_ent_tri = libObject.torus_1d_entropy_triangular
    torus_1d_ent_tri.restype = ctypes.c_double
    torus_1d_ent_tri.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return torus_1d_ent_tri(n, r0, lim)

def torus_1d_entropy_exclusion(n, r0, lim):
    torus_1d_ent_excl = libObject.torus_1d_entropy_exclusion
    torus_1d_ent_excl.restype = ctypes.c_double
    torus_1d_ent_excl.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return torus_1d_ent_excl(n, r0, lim)


def torus_2d_entropy_rayleigh(n, r0, eta, lim):
    torus_2d_ent_ray = libObject.torus_2d_entropy_rayleigh
    torus_2d_ent_ray.restype = ctypes.c_double
    torus_2d_ent_ray.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]

    return torus_2d_ent_ray(n, r0, eta, lim)

def torus_2d_entropy_hard(n, r0, lim):
    torus_2d_ent_hard = libObject.torus_2d_entropy_hard
    torus_2d_ent_hard.restype = ctypes.c_double
    torus_2d_ent_hard.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return torus_2d_ent_hard(n, r0, lim)

def torus_2d_entropy_triangular(n, r0, lim):
    torus_2d_ent_tri = libObject.torus_2d_entropy_triangular
    torus_2d_ent_tri.restype = ctypes.c_double
    torus_2d_ent_tri.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return torus_2d_ent_tri(n, r0, lim)

def torus_2d_entropy_exclusion(n, r0, lim):
    torus_2d_ent_excl = libObject.torus_2d_entropy_exclusion
    torus_2d_ent_excl.restype = ctypes.c_double
    torus_2d_ent_excl.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return torus_2d_ent_excl(n, r0, lim)

def line_entropy_rayleigh(n, r0, eta, lim):
    line_ent_ray = libObject.line_entropy_rayleigh
    line_ent_ray.restype = ctypes.c_double
    line_ent_ray.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]

    return line_ent_ray(n, r0, eta, lim)

def line_entropy_hard(n, r0, lim):
    line_ent_hard = libObject.line_entropy_hard
    line_ent_hard.restype = ctypes.c_double
    line_ent_hard.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return line_ent_hard(n, r0, lim)

def line_entropy_triangular(n, r0, lim):
    line_ent_tri = libObject.line_entropy_triangular
    line_ent_tri.restype = ctypes.c_double
    line_ent_tri.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return line_ent_tri(n, r0, lim)

def line_entropy_exclusion(n, r0, lim):
    line_ent_excl = libObject.line_entropy_exclusion
    line_ent_excl.restype = ctypes.c_double
    line_ent_excl.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return line_ent_excl(n, r0, lim)


def disc_entropy_rayleigh(n, r0, eta, lim):
    disc_ent_ray = libObject.disc_entropy_rayleigh
    disc_ent_ray.restype = ctypes.c_double
    disc_ent_ray.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]

    return disc_ent_ray(n, r0, eta, lim)

def disc_entropy_hard(n, r0, lim):
    disc_ent_hard = libObject.disc_entropy_hard
    disc_ent_hard.restype = ctypes.c_double
    disc_ent_hard.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return disc_ent_hard(n, r0, lim)

def disc_entropy_triangular(n, r0, lim):
    disc_ent_tri = libObject.disc_entropy_triangular
    disc_ent_tri.restype = ctypes.c_double
    disc_ent_tri.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return disc_ent_tri(n, r0, lim)

def disc_entropy_exclusion(n, r0, lim):
    disc_ent_excl = libObject.disc_entropy_exclusion
    disc_ent_excl.restype = ctypes.c_double
    disc_ent_excl.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int]

    return disc_ent_excl(n, r0, lim)



# p_bar integrators

def p_bar_hard(geom, r0, lim):
    geom = geom.encode('utf-8')
    p_bar_hard = libObject.p_bar_hard
    p_bar_hard.restype = ctypes.c_double
    p_bar_hard.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_int]

    return p_bar_hard(geom, r0, lim)

def p_bar_triangular(geom, r0, lim):
    geom = geom.encode('utf-8')
    p_bar_tri = libObject.p_bar_triangular
    p_bar_tri.restype = ctypes.c_double
    p_bar_tri.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_int]

    return p_bar_tri(geom, r0, lim)

def p_bar_rayleigh(geom, r0, eta, lim):
    geom = geom.encode('utf-8')
    p_bar_ray = libObject.p_bar_rayleigh
    p_bar_ray.restype = ctypes.c_double
    p_bar_ray.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_double, ctypes.c_int]

    return p_bar_ray(geom, r0, eta, lim)

def p_bar_exclusion(geom, r0, lim):
    geom = geom.encode('utf-8')
    p_bar_excl = libObject.p_bar_exclusion
    p_bar_excl.restype = ctypes.c_double
    p_bar_excl.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_int]

    return p_bar_excl(geom, r0, lim)


def entropy_sim(n, geom, cf, r0, eta=2, lim=100000):

    if geom == "square":
        if cf == "hard":
            return square_entropy_hard(n, r0, lim)
        elif cf == "triangular":
            return square_entropy_triangular(n, r0, lim)
        elif cf == "rayleigh":
            return square_entropy_rayleigh(n, r0, eta, lim)
        elif cf == "waxman":
            return square_entropy_rayleigh(n, r0, 1, lim)
        elif cf == "exclusion":
            return square_entropy_exclusion(n, r0, lim)
        else:
            print("Invalid correlation function, please choose from 'hard', 'triangular', 'rayleigh'")
    elif geom == "torus_1d":
        if cf == "hard":
            return torus_1d_entropy_hard(n, r0, lim)
        elif cf == "triangular":
            return torus_1d_entropy_triangular(n, r0, lim)
        elif cf == "rayleigh":
            return torus_1d_entropy_rayleigh(n, r0, eta, lim)
        elif cf == "waxman":
            return torus_1d_entropy_rayleigh(n, r0, 1, lim)
        elif cf == "exclusion":
            return torus_1d_entropy_exclusion(n, r0, lim)
    elif geom == "torus_2d":
        if cf == "hard":
            return torus_2d_entropy_hard(n, r0, lim)
        elif cf == "triangular":
            return torus_2d_entropy_triangular(n, r0, lim)
        elif cf == "rayleigh":
            return torus_2d_entropy_rayleigh(n, r0, eta, lim)
        elif cf == "waxman":
            return torus_2d_entropy_rayleigh(n, r0, 1, lim)
        elif cf == "exclusion":
            return torus_2d_entropy_exclusion(n, r0, lim)
    elif geom == "line":
        if cf == "hard":
            return line_entropy_hard(n, r0, lim)
        elif cf == "triangular":
            return line_entropy_triangular(n, r0, lim)
        elif cf == "rayleigh":
            return line_entropy_rayleigh(n, r0, eta, lim)
        elif cf == "waxman":
            return line_entropy_rayleigh(n, r0, 1, lim)
        elif cf == "exclusion":
            return line_entropy_exclusion(n, r0, lim)
    elif geom == "disc":
        if cf == "hard":
            return disc_entropy_hard(n, r0, lim)
        elif cf == "triangular":
            return disc_entropy_triangular(n, r0, lim)
        elif cf == "rayleigh":
            return disc_entropy_rayleigh(n, r0, eta, lim)
        elif cf == "waxman":
            return disc_entropy_rayleigh(n, r0, 1, lim)
        elif cf == "exclusion":
            return disc_entropy_exclusion(n, r0, lim)
    else:
        print("Invalid geometry, please choose from 'square'")


def f_square(r0):
    f_square = libObject.f_square
    f_square.restype = ctypes.c_double
    f_square.argtypes = [ctypes.c_double]

    return f_square(r0)
