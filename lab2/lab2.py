import math
import numpy as np
import pylab
import scipy.integrate
from decimal import Decimal

# extended trapezoidal integration approximation of func 
# from a to b in N subintervals
def extrap(func, a, b, N):	
    a = Decimal(a)
    b = Decimal(b)
    N = Decimal(N)
    hN = (b - a)/N # width of subinterval
    x = np.arange(a, (N + 1) * hN, hN) # array from a to b separated by hN
    fx = func(x)
    
    return (np.sum(fx) - fx[0]/Decimal(2) - fx[-1]/Decimal(2)) * hN

# extended Simpson's integration approximation of func
# from a to b in N subintervals
def exsimp(func, a, b, N):
    a = Decimal(a)
    b = Decimal(b)
    N = Decimal(N)
    hN = (b - a)/N # width of subinterval 
    x = np.arange(a, b + hN/2, hN/2) # array from a to b separated by hN/2
    fx = func(x)

    fa = fx[0::2] # f(a and b values)
    fc = fx[1::2] # f(middle elements)
    
    return (np.sum(fa)/Decimal(3) - fa[0]/Decimal(6) - fa[-1]/Decimal(6) + np.sum(fc)*Decimal(4.0/6))*hN

# refine Simpson's for N = N0, 2N0,... 
# until |I_simp(2^k*N0)-I_simp(2^(k+1)*N|/I_simp(2^kN0) < error
def refinesimp(func, a, b, e):
    N0 = 1
    i = 0 # index
    error = 1
    
    while error >= e:
        S1 = exsimp(func, a, b, (2**i)*N0)
        i += 1
        S2 = exsimp(func, a, b, (2**(i))*N0)
        error = math.fabs(S1-S2)/S1

        print S1
        print S2
        print error

    return S1

N = 2 * np.rint(np.logspace(1.0, 3.0, num=100))
Itrap = np.array([])
Isimp = np.array([])
I = Decimal(np.e - 1) # true value of integral

for i in N:
    Itrap = np.append(Itrap, extrap(np.exp, 0, 1, i))
    Isimp = np.append(Isimp, exsimp(np.exp, 0, 1, i))

#pylab.scatter(N, Itrap-I)
#pylab.scatter(N, Isimp-I)
#pylab.xlabel('N')
#pylab.ylabel('I_appr-I')
#pylab.show()

#print Itrap-I

pylab.loglog(N, Itrap-I)
pylab.loglog(N, Isimp-I)
pylab.xlabel('N')
pylab.ylabel('I_appr - I')
pylab.show()

'''
print extrap(np.exp, 0, 1, 100) - I # 1.43189913717e-05
print exsimp(np.exp, 0, 1, 100) - I # 5.96589444513e-12

Is = refinesimp(np.exp, 0, 1, 0.000001)
print Is # 1.71831884192
print (Is - I) # 3.70134627019e-05

# (1.7182818284590453, 1.9076760487502457e-14)
print scipy.integrate.quad(np.exp, 0, 1) 
# 1.71828182846
print scipy.integrate.romberg(np.exp, 0, 1)
'''
def square (x):
    return (x*x)

#I1 = 0.333333333333 # integral of x^2 from 1 to 10
#Is1 = refinesimp(square, 0, 1, 0.000001)

#print Is1 # 0.333333333333
#print (Is1 - I1) # 1.66533453694e-16
