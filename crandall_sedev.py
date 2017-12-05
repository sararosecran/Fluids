#1D Sedev blast wave

import numpy as np
import matplotlib.pyplot as plt


def velocity(steps, dr, rho0, P, q, n, u):
	u_new = np.zeros(steps)
	#courant condition
	m = np.asarray(steps)
	for i in range(2, steps-2):
		m = np.append(m, abs(u[n-1][i]) + (gamma * P[n-1][i] * v[n-1][i]) / dr)
	dt = 0.5 * np.min(m)
	print dt
	for j in range(2, steps - 2):
		u_new[j-1] = - dt / rho0 * (P[n-1][j] - P[n-1][j-2] + q[n-1][j] - q[n-1][j-2]) / dr + u[n-2][j-1]
	return u_new
    
def radius(steps, r_array, dt, u, n):
	r_new = np.zeros(steps)
	for j in range(2, steps - 2):
	    r_new[j-1] = u[n][j-1] * dt + r_array[n-1][j-1]
	return r_new

def inverse_rho(steps, rho0, r_array, dr):
	v_new = np.zeros(steps)
	for j in range(2, steps-2):
		v_new[j] = 1.0 / rho0 * (r_array[n+1][j+1] - r_array[n+1][j-1]) / dr
	return v_new

def art_vis(steps, u, v, a):
	q_new = np.zeros(steps)
	for j in range(2, steps-2):
		if (u[n][j+1] - u[n][j-1]) < 0. :
			q_new[j] = (0.5 * a**2 * (u[n][j+1]-u[n][j-1])**2) / (v[n+1][j] + v[n-1][j])
		else:
			q_new[j] = 0
	return q_new

def energy(steps, P, q, v, E):
	e_new = np.zeros(steps)
	for j in range(2, steps-2):
		e_new[j] = -(P[n-1][j] + q[n][j]) * (v[n+1][j] - v[n-1][j]) + E[n-1][j]
	return e_new

def pressure(steps, E, v, gamma):
	P_new = np.zeros(steps)
	for j in range(2, steps-2):
		P_new[j] = E[n+1][j] * (gamma - 1.0) / v[n+1][j]
	return P_new

dr = 1000.0
dt = 1000.0
P0 = 0.7 # bars = 70KPa
gamma = 5.0 / 3.0 #monotomic gas
rho0 = 0.001 #g/cm^3
a = (gamma * P0 / rho0)**0.5
steps = 10
u = np.zeros((steps,steps))
u[0] = u[1] = np.full(steps,0.0) #inital velocity is zero
P = np.zeros((steps,steps))
for i in range(steps): #inital pressure. All r<8 will have a large Pressure.
	if i < 3:
		P[0][i] = P[1][i] = P0
	else:
		P[0][i] = P[1][i] = 0.1
q = np.zeros((steps, steps))
q[0] = q[1] = np.full(steps,0.0) #inital viscosity is zero
r_array = np.zeros((steps,steps))
v = np.zeros((steps, steps))
E = np.zeros((steps, steps))
E[0] = E[1] = P0 * np.pi * (3.0 *dr)**2 # initial Energy inside blast




#step in time
for n in range(2,steps-2):
	u[n] = velocity(steps, dr, rho0, P, q, n, u)
	r_array[n+1] = radius(steps, r_array, dt, u, n)
	v[n+1] = inverse_rho(steps, rho0, r_array, dr)
	q[n] = art_vis(steps, u, v, a)
	E[n+1] = energy(steps, P, q, v, E )
	P[n+1] = pressure(steps, E, v, gamma)

print u
#plt.plot(steps, u[0])
#plt.show()