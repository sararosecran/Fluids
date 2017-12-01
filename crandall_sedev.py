#1D Sedev blast wave


def velocity(r_arraydx,dt,rho0,P,q):
    u_new = np.zeros(len(r_array))
	for j in range(2, len(r_array - 1):
        u_new[j] = - dt / rho0 * (P[j] - P[j-2] + q[j] - q[j-2]) / dx + u[j-1]
	return u_new
    
def x(r_array):
	x_new = np.zeros(len(r_array))
    for j in range(2, len(r_array - 1))
	    


#for i in range(

dr = 1.0
dt = 1.0
P0 = 0.7 # bars = 70KPa
rho0 = 0.001 #g/cm^3
r_array = np.arange(1.0,100.0,dr)
