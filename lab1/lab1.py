#!/usr/bin/python

import math, sys
import numpy as np
import pylab

fx = float(sys.argv[1]) # x frequency
fy = float(sys.argv[2]) # y frequency
Ax = float(sys.argv[3]) # x amplitude
Ay = float(sys.argv[4]) # y amplitude
phi = float(sys.argv[5]) # phase difference
dt = float(sys.argv[6]) # step size
N = int(sys.argv[7]) # number of steps


[x, y, z] = [[], [], []] # using lists in python
t = 0

while t <= N * dt:
    newx = Ax * math.cos(2 * math.pi * fx * t)
    newy = Ay * math.sin(2 * math.pi * fy * t + phi)
    x.append(newx)
    y.append(newy)
    z.append(newx + newy)
    t += dt

file = open('Output.txt', 'w')

file.write('X\tY\tZ\n')

# writes x, y, z into Output.txt with tabs separating columns
for i in range(len(x)):
    file.write(str(x[i]) + '\t' + str(y[i]) + '\t' + str(z[i]) + '\n')

file.close()

# using arrays in numpy
 
t = np.arange(0, (N + 1) * dt, dt)

x = Ax * np.cos(2 * np.pi * fx * t) # operates on each element of array
y = Ay * np.sin(2 * np.pi * fy * t + phi)
z = x + y
'''
# writes x, y, z into three different files, one for each variable
np.savetxt('X_Output.txt', x, delimiter='\n')
np.savetxt('Y_Output.txt', y, delimiter='\n')
np.savetxt('Z_Output.txt', z, delimiter='\n')
'''
plot = raw_input('Plot (XY/ZT): ')

if (plot == 'XY'):
    pylab.plot(x, y)
    pylab.xlabel('X')
    pylab.ylabel('Y')

elif (plot == 'ZT'):
    pylab.plot(t, z)
    pylab.xlabel('t')
    pylab.ylabel('Z(t)')

pylab.show()
