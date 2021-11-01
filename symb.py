from sympy import *
import mpmath
from mpmath import mp
mp.dps=100
import numpy as np
k1,k2,a,t=symbols("k1,k2,a,t",cls=Symbol)
#
# k1=0.1
# k2=0.2
k1s=conjugate(k1)
k2s=conjugate(k2)

ebeta3=1/(exp(k2+k1s)+exp(-k2-k1s)-2)
ebeta1=1/(exp(k1+k1s)+exp(-k1-k1s)-2)
ebeta2=1/(exp(k1+k2s)+exp(-k1-k2s)-2)
ebeta4=1/(exp(k2+k2s)+exp(-k2-k2s)-2)



##the answer
rho1=Rational(1,4)/(sinh((k1+k1s)/2))**2
rho3=Rational(1,4)/(sinh((k2+k1s)/2))**2
chi=exp(k1-k2)+exp(k2-k1)-2#4*(sinh((k1-k2)/2))**2
chis=conjugate(chi)
rho2=Rational(1,4)*1/(sinh((k1+k2s)/2))**2
rho4=Rational(1,4)*1/(sinh((k2+k2s)/2))**2
tau1=rho1*rho3*chi
tau2=rho2*rho4*chi

d1=exp(k1+k1s)+exp(-k1-k1s)-2
d2=exp(k1+k2s)+exp(-k1-k2s)-2
d3=exp(k2+k1s)+exp(-k2-k1s)-2
d4=exp(k2+k2s)+exp(-k2-k2s)-2

dm1=-exp(k1-I*a*t)+exp(-k1+I*a*t)
dm2=-exp(k2-I*a*t)+exp(-k2+I*a*t)
dm1s=-exp(k1s+I*a*t)+exp(-k1s-I*a*t)
dm2s=-exp(k2s+I*a*t)+exp(-k2s-I*a*t)



N1x=-chis*(dm2+dm1s+dm2s)+d4*(dm2+dm1s)\
    +d3*(dm2+dm2s)-d4*dm2s-d3*dm1s\
    +chis*exp(-I*a*t-k2-k1s-k2s)+d4*exp(-I*a*t+k2+k1s-k2s)\
    +d3*exp(-I*a*t+k2+k2s-k1s)-chis*exp(I*a*t+k2+k1s+k2s)\
    -d4*exp(I*a*t-k2-k1s+k2s)-d3*exp(I*a*t-k2+k1s-k2s)

N2x=-chis*(dm1+dm1s+dm2s)+d2*(dm1+dm1s)\
    +d1*(dm1+dm2s)-d2*dm2s-d1*dm1s\
    +chis*exp(-I*a*t-k1-k1s-k2s)+d2*exp(-I*a*t+k1+k1s-k2s)\
    +d1*exp(-I*a*t+k1+k2s-k1s)-chis*exp(I*a*t+k1+k1s+k2s)\
    -d2*exp(I*a*t-k1-k1s+k2s)-d1*exp(I*a*t-k1+k1s-k2s)

rst=a/t+t/a
