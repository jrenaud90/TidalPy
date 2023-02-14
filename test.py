from numba import njit
from scipy import special

def f():
    print('Gamma', special.gamma(1.2))

fn = njit(f)

f()
fn()