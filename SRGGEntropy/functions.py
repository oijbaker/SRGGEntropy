import numpy as np

def H2(p):
    return -p*np.log(p) - (1-p)*np.log(1-p)