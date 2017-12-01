#1D Sedev blast wave

import numpy as np


def velocity(r_array,dr,dt,rho0,P,q,n,u):
	u_new = np.zeros(len(r_array))
	for j in range(2, len(r_array - 2)):
		u_new[j-1] = - dt / rho0 * (P[n-1][j] - P[n-1][j-2] + q[n-1][j] - q[n-1][j-2]) / dr #+ u[n-2][j-1]
	return u_new
    
def radius(r_array,dt,u,r,n):
	r_new = np.zeros(len(r_array))
	for j in range(2, len(r_array - 2)):
	    r_new[j-1] = u[n][j-1] * dt + r[n-1][j-1]
	return r_new



dr = 1.0
dt = 1.0
P0 = 0.7 # bars = 70KPa
rho0 = 0.001 #g/cm^3
r_array = np.arange(1.0,101.0,dr)
u = np.zeros([len(r_array),len(r_array)])
u[0] = u[1] = np.full(len(r_array),0.0) #inital velocity is zero
P = np.zeros([len(r_array),len(r_array)])
for i in range(len(r_array)): #inital pressure. All r<8 will have a large Pressure.
	if i < 8:
		P[0][i] = P[1][i] = P0
	else:
		P[0][i] = P[1][i] = 0.1
q = np.zeros([len(r_array),len(r_array)])
q[0] = q[1] = np.full(len(r_array),0.0) #inital viscosity is zero


#step in time
for n in range(2,len(r_array)-2):
	u[n] = velocity(r_array, dr, dt, rho0, P, q, n, u)
	r[n+1] = radius(r_array,
