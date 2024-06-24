import ctypes

libObject = ctypes.CDLL('C:\\Users\\dh23880\\Documents\\SRGGEntropy-1\\SRGGEntropy\\clibrary.so')

sq_ent_hard = libObject.square_entropy_hard
sq_ent_hard.restype = ctypes.c_double
sq_ent_hard.argtypes = [ctypes.c_int, ctypes.c_double,ctypes.c_int]

sq_ent_ray = libObject.square_entropy_rayleigh
sq_ent_ray.restype = ctypes.c_double
sq_ent_ray.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]

print(sq_ent_ray(3, 0.61, 2, 1000000))