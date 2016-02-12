import numpy as np
import pylab

def exspring_anal(x0, v0, tf, h, k, m):
    # uses explicit Euler method to determine
    # (positions, velocities, times) tuple of spring
    # spring constant = k, mass = m
    # initial position = x0, initial velocity = v0, initial time t0 = 0
    # final time evaluated at tf
    # step size = h

    N = tf/h # total number of steps required
    t = np.arange(0, tf+h, h) # array of times

    x = np.array([x0]) # begin with the initial position, velocity
    v = np.array([v0])

    for ti in t[1:]:
        xi = x[-1] + h * v[-1]
        vi = v[-1] - h * (k / float(m)) * x[-1]
	
        x = np.append(x, xi)
        v = np.append(v, vi)

    return (x, v, t)

def imspring_anal(x0, v0, tf, h, k, m):
    # uses implicit Euler method to determine
    # (positions, velocities, times) tuple of spring
    # spring constant = k, mass = m
    # initial position = x0, initial velocity = v0, initial time t0 = 0
    # final time evaluated at tf
    # step size = h

    N = tf/h
    t = np.arange(0, tf+h, h)

    x = np.array([x0])
    v = np.array([v0])

    for ti in t[1:]:
        xi = (x[-1] + h * v[-1])/(float(h * h) + 1) # determines next position
        vi = v[-1] - h * (k / float(m)) * xi

        x = np.append(x, xi)
        v = np.append(v, vi)

    return (x, v, t)

def spring_xerror(x0, v0, tf, h, k, m):
    # finds maximum error in analytic approximation of position
    
    (xi, vi, t) = exspring_anal(x0, v0, tf, h, k, m)

    # exact solutions of spring
    w = np.sqrt(float(k)/m)
    x = np.sin(w * t)
    v = w * np.cos(w * t)

    return np.amax(np.abs(x - xi))

# simulate spring with k = m = 1 for t = 0 to 100
(xi, vi, t) = exspring_anal(0, 1, 100, 0.01, 1, 1)
E_e = xi**2 + vi**2

(ixi, ivi, t) = imspring_anal(0, 1, 100, 0.01, 1, 1)
E_i = ixi**2 + ivi**2

#pylab.plot(t, xi)
pylab.plot(t, ixi)
pylab.xlabel('Time')
pylab.ylabel('Position')
pylab.show()

#pylab.plot(t, vi)
pylab.plot(t, ivi)
pylab.xlabel('Time')
pylab.ylabel('Velocity')
pylab.show()

# exact solutions of spring for k = m = 1 with IC x0 = 0, v0 = 1
x = np.sin(t)
v = np.cos(t)

e_x = np.abs(x - xi)
e_v = np.abs(v - vi)

e_ix = np.abs(x - ixi)
e_iv = np.abs(v - ivi)

#pylab.plot(t, e_x)
#pylab.plot(t, np.abs(E_e - 1))
#pylab.plot(t, e_v)
#pylab.plot(t, e_ix)
#pylab.plot(t, e_iv)
#pylab.plot(t, np.abs(E_i - 1))
#pylab.xlabel('Time')
#pylab.ylabel('Error')
#pylab.show()

# finds error for different step sizes (x0 = 0, v0 = 1, tf = 10, k = m= 1)
h = 0.1/np.power(2, np.arange(0, 10))
e = np.array([])

for hi in h:
    e = np.append(e, spring_xerror(0, 1, 10, hi, 1, 1))

#pylab.plot(h, e)
#pylab.xlabel('Step Size (h)')
#pylab.ylabel('Error')
#pylab.show()

#pylab.plot(t, E_e)
pylab.plot(t, E_i)
pylab.xlabel('Time')
pylab.ylabel('Energy')
pylab.show()
