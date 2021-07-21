from scipy.misc import derivative
import numpy as np

from sympy import *

bt=0.4
N=2**10
a=16/N
x0=N/2

n0=N/2

def x(t):
    return 2/(a*bt)*np.sinh(bt)*(np.cos(a*t)-1)+x0

def phi(t):
    return -2/a*np.cosh(bt)*np.sin(a*t)

def psi(n,t):
    return np.sinh(bt)*1/np.cosh(bt*(n-x(t)))*np.exp(-1j*(phi(t)+a*n*t))

def psiN(t):
    return psi(n0,t)

def psiNp1(t):
    return psi(n0+1,t)

def psiNm1(t):
    return psi(n0-1,t)



t0=10

lhs=1j*derivative(psiN,t0,dx=1e-6)

rhs=-(psiNp1(t0)+psiNm1(t0))*(1+np.abs(psiN(t0))**2)+a*n0*psiN(t0)

diff=lhs-rhs
print(lhs)
print(diff)
