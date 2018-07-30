import numpy as np


def D(T,material):
    R=8.314
    k_b=8.6e-5
    if material=="steel":
        return 7.3e-7*np.exp(-6.3e3/T)
    elif material=="polymer":
        return 2.0e-7*np.exp(-29000.0/R/T)
    elif material=="concrete":
        return 1e-6
    elif material=="tungsten":
        return 4.1e-7*np.exp(-0.39/k_b/T)
    elif material=="lithium_lead":
        return 2.5e-7*np.exp(-27000.0/R/T)

def heat_capacity(T,material):
    return
def thermal_conductivity(T,material):
    return
def specific_heat(T,material):
    return
def density(T,material):
    return
