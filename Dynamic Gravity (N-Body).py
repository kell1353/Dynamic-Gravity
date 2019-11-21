from mayavi import mlab
from matplotlib import pyplot as plt
import numpy as np
import math as m
import time
import copy

mlab_fig = mlab.figure('N - Body System', bgcolor = (0,0,0), size = (700,500))

'Spherical constants'
points_range = 50
phi = np.linspace(0, 2*np.pi, points_range)
theta = np.linspace(0, 2*np.pi, points_range)
phi, theta = np.meshgrid(phi, theta)

'Common Variables'
sqrt = m.sqrt
sin, cos, tan = np.sin, np.cos, np.tan
pi = np.pi

def draw_object(r, x, y, z):
    x = (r* np.sin(phi) * np.cos(theta)) + x
    y = (r*np.sin(phi) * np.sin(theta)) + y
    z = (r*np.cos(phi)) + z
    sphere = mlab.mesh(x, y, z)

def draw_sphere(e_r, p_r):
    global x; global y; global z
    x = (e_r* np.sin(phi) * np.cos(theta))
    y = (e_r*np.sin(phi) * np.sin(theta))
    z = (p_r*np.cos(phi))

'Calculate the distance of a vector to the closest point on a line'
def vectorMag(vx, vy, vz):
    return sqrt(vx**2 + vy**2 + vz**2)

def vec(x0, y0, z0, x1, y1, z1):
    global ux, uy, uz
    'Calculate all of the Relevent Vectors'
    ux, uy, uz = x1 - x0, y1 - y0, z1 - z0

def dist(x0, y0, z0, x1, y1, z1):
    return sqrt((x0 - x1)**2 + (y0 - y1)**2 + (z0 - z1)**2)

'Draw a vector fields given the tail locations and components'
def draw_velVector(x0, y0, z0, x1, y1, z1):
    global vx, vy, vz
    vec(x0, y0, z0, x1, y1, z1)
    vx, vy, vz = ux, uy, uz
    mlab.quiver3d(x0, y0, z0, vx, vy, vz, line_width = 1, scale_factor= 1)

def draw_gravVector(M, x1, y1, z1):
    global fx, fy, fz
    x0, y0, z0 = 0, 0, 0
    r = dist(x0, y0, z0, x1, y1, z1)
    G = 1
    f_g = (G*M)/r
    vec(x1, y1, z1, x0, y0, z0)
    vecNorm = vectorMag(ux, uy, uz)
    fx, fy, fz = f_g*(ux/vecNorm), f_g*(uy/vecNorm), f_g*(uz/vecNorm)
    mlab.quiver3d(x1, y1, z1, fx, fy, fz, color = (0, 1, 0), line_width = 1, scale_factor= .000000001)

#https://sites.temple.edu/math5061/files/2016/12/final_project.pdf
'x0, y0, z0 are the cooridnates of object with the gravity'
'x1, y1, z1 are the cooridnates of object your trying to calculate'
'M is the mass of the large object'
def a(M, x0, y0, z0, x1, y1, z1):
    global ax; global ay; global az    
    del_x, del_y, del_z = x0 - x1, y0 - y1, z0 - z1
    d = del_x**2 + del_y**2 + del_z**2
##    Fx = ((G*M*m)/(d))*(del_x/(sqrt(d)))
##    Fy = ((G*M*m)/(d))*(del_y/(sqrt(d)))
##    Fz = ((G*M*m)/(d))*(del_z/(sqrt(d)))
    
    ax = ((G*M)/(d))*(del_x/(sqrt(d)))
    ay = ((G*M)/(d))*(del_y/(sqrt(d)))
    az = ((G*M)/(d))*(del_z/(sqrt(d)))
    #print(ax, ay, az)

def a_sum(ID, xp, yp, zp):
    global ax_sum; global ay_sum; global az_sum  
    ax_sum, ay_sum, az_sum = 0, 0, 0
    for i in range(0, len(u)):
        if i == ID:
            ax_sum += 0; ay_sum += 0; az_sum += 0
        else:
            if instance == 'animate':
                a(v[i][1], v[i][4], v[i][5], v[i][6], xp, yp, zp)
                ax_sum += ax; ay_sum += ay; az_sum += az
            if instance == 'trajectories':
                a(u[i][1], u[i][4], u[i][5], u[i][6], xp, yp, zp)
                ax_sum += ax; ay_sum += ay; az_sum += az

'Compute the next x,y,z positions and velocities using the 4th-order Runge-Kutta technique'
def RK_4(ID, x, y, z, vx, vy, vz):
    global x_n1; global y_n1; global z_n1; global vx_n1; global vy_n1; global vz_n1
    'Step Length: 1 day = 86,400 (in seconds (s))'
    del_t = 86400

    k, kv = {}, {}
    for i in range(1,5):
        if i == 1:
            k['k1_x'] = vx
            k['k1_y'] = vy
            k['k1_z'] = vz
            a_sum(ID, x, y, z)
        if i == 2 or i == 3:
            k['k'+str(i)+'_x'] = vx + ((del_t/2)*kv['k'+str(i-1)+'_vx'])
            k['k'+str(i)+'_y'] = vy + ((del_t/2)*kv['k'+str(i-1)+'_vy'])
            k['k'+str(i)+'_z'] = vz + ((del_t/2)*kv['k'+str(i-1)+'_vz'])
            a_sum(ID, x + ((del_t/2)*k['k'+str(i-1)+'_x']), y + ((del_t/2)*k['k'+str(i-1)+'_y']), z + ((del_t/2)*k['k'+str(i-1)+'_z']))
        if i == 4:
            k['k'+str(i)+'_x'] = vx + ((del_t)*kv['k'+str(i-1)+'_vx'])
            k['k'+str(i)+'_y'] = vy + ((del_t)*kv['k'+str(i-1)+'_vy'])
            k['k'+str(i)+'_z'] = vz + ((del_t)*kv['k'+str(i-1)+'_vz'])
            a_sum(ID, x + ((del_t)*k['k3_x']), y + ((del_t)*k['k3_y']), z + ((del_t)*k['k3_z']))
            
        kv['k'+str(i)+'_vx'] = ax_sum
        kv['k'+str(i)+'_vy'] = ay_sum
        kv['k'+str(i)+'_vz'] = az_sum


    x_n1 = x + ((del_t/6)*(k['k1_x'] + (2*k['k2_x']) + (2*k['k3_x']) + k['k4_x']))
    y_n1 = y + ((del_t/6)*(k['k1_y'] + (2*k['k2_y']) + (2*k['k3_y']) + k['k4_y']))
    z_n1 = z + ((del_t/6)*(k['k1_z'] + (2*k['k2_z']) + (2*k['k3_z']) + k['k4_z']))

    vx_n1 = vx + ((del_t/6)*(kv['k1_vx'] + (2*kv['k2_vx']) + (2*kv['k3_vx']) + kv['k4_vx']))
    vy_n1 = vy + ((del_t/6)*(kv['k1_vy'] + (2*kv['k2_vy']) + (2*kv['k3_vy']) + kv['k4_vy']))
    vz_n1 = vz + ((del_t/6)*(kv['k1_vz'] + (2*kv['k2_vz']) + (2*kv['k3_vz']) + kv['k4_vz']))

    

'Set Initial Positions (in meters (m)) and Set Initial Velocities (in m/s) and Mass (in kilograms (kg)) of the gravitational objects'
'gravObjects = [(ID, Mass, EqRadius, PolarRadius, Initial x_pos, Initial y_pos, Initial z_pos, Initial x_vel, Initial y_vel, Initial z_vel)]'
#for i in range(7750, 7851, 10):
mlab_fig = mlab.figure('N - Body System', bgcolor = (0,0,0), size = (700,500))
#print(i)

#Sun, Earth and Mars
##u = [[0, 1.989*(10**30), 6955100000, 6955100000, 0, 0, 0, 0, 0, 0], \
##         [1, 5.972*(10**24), 237100000, 237100000, 152100000000, 0, 0, 0, 29290, 0], \
##         [2, 7.348*(10**22), 173710000, 173710000, 152484400000, 0, 0, 0, 29451, 0]]
##u = [[0, 1.989*(10**30), 16955100000, 16955100000, -20000000000*3, 0, 0, 0, 12000, 0, (1, 0, 0)], \
##         [1, 1.989*(10**30), 16955100000, 16955100000, 115000000000*3, 78000000000*2, 0, -17000, -5000, 0, (0, 0, 1)], \
##         [2, 1.989*(10**30), 16955100000, 16955100000, 115000000000*3, -78000000000*2, 0, 12000, -1000, 0, (0, 1, 0)]]


'Three Body Solution (Figure Eight)'
V = 8850
v0x, v0y  = -2*(V*.93420737), -2*(V*.864731) 
v1y, v1y = V*.93420737, V*.864731

u = [[0, 1.989*(10**30), 16955100000, 16955100000, 0, 0, 0, v0x, v0y, 0, (1, 0, 0)], \
         [1, 1.989*(10**30), 16955100000, 16955100000, 545000000000*.97000436, -545000000000*.24308753, 0, v1y, v1y, 0, (0, 0, 1)], \
         [2, 1.989*(10**30), 16955100000, 16955100000, -545000000000*.97000436, 545000000000*.24308753, 0, v1y, v1y, 0, (0, 1, 0)]]

v = copy.deepcopy(u)

'Set Initial Conditions'
'Gravitational Constant (in N*m^2/(kg^2))'
G = 6.67408*(10**(-11))

#FOR TESTING THE RK4 FUNCTION
##instance = 'trajectories'
##u = [[0, 1.989*(10**30), 16955100000, 16955100000, 0, 0, 0, 0, 00, 0, (1, 0, 1)], \
##         [1, 5.972*(10**24), 1955100000, 1955100000, -147095000000, 0, 0, 0, -30300, 0, (0, 0, 1)]]
##RK_4(u[1][0], u[1][4], u[1][5], u[1][6], u[1][7], u[1][8], u[1][9])
##print(x_n1, y_n1, z_n1)
##print(vx_n1, vy_n1, vz_n1)


'Define Initial Objects using values from 'u' and Create Dictionary to house them'
d_obj = {}
for i in range(len(u)):
    draw_sphere(u[i][2], u[i][3])
    d_obj['object'+str(i)] = mlab.mesh(x, y, z, color = u[i][10])

'Define Initial Trajectories and Create Dictionary to house them'
d_traj = {}
objDist = {}
for i in range(len(u)):
    d_traj['x_traj'+str(i)] = []
    d_traj['y_traj'+str(i)] = []
    d_traj['z_traj'+str(i)] = []
    objDist['objectDist'+str(i)] = []

days = 1000
'Draw the Trajectory Lines for each objects updating for each objects position'
instance = 'trajectories'
index = []
#objectDist, object2Dist, object3Dist, , [], [], []
for i in range(0, days):
    index.append(i)
    for j in range(0, len(u)):
        RK_4(u[j][0], u[j][4], u[j][5], u[j][6], u[j][7], u[j][8], u[j][9])
        d_traj['x_traj'+str(j)].append(x_n1)
        d_traj['y_traj'+str(j)].append(y_n1)
        d_traj['z_traj'+str(j)].append(z_n1)
        objDist['objectDist'+str(j)].append(dist(x_n1, y_n1, z_n1, 0, 0, 0))
        u[j][4], u[j][5], u[j][6] = x_n1, y_n1, z_n1
        u[j][7], u[j][8], u[j][9] = vx_n1, vy_n1, vz_n1

for i in range(len(u)):
    mlab.plot3d(d_traj['x_traj'+str(i)], d_traj['y_traj'+str(i)], d_traj['z_traj'+str(i)], color=(1,1,1), tube_radius = 100000000/2)


'Animate the figure'
instance = 'animate'
@mlab.animate(delay = 10)
def anim():
    for i in range(0, days):
        #print(i)
        time.sleep(.00000000000000001)
        for j in range(0, len(v)):
            RK_4(v[j][0], v[j][4], v[j][5], v[j][6], v[j][7], v[j][8], v[j][9])            
            d_obj['object'+str(j)].mlab_source.set(x = v[j][2] * np.sin(phi) * np.cos(theta) + x_n1, y = v[j][2] * np.sin(phi) * np.sin(theta) + y_n1, z = v[j][3] * np.cos(phi) + z_n1)
                
            v[j][4], v[j][5], v[j][6] = x_n1, y_n1, z_n1
            v[j][7], v[j][8], v[j][9] = vx_n1, vy_n1, vz_n1
        yield
    print(i)
anim()
    
mlab.view(distance=1250000000000, focalpoint = (0, 0, 0))
mlab.show()
#mlab.draw()


fig = plt.figure()
ax = fig.add_subplot(211)
ax.set_xlabel('Days')
ax.set_ylabel('Distance from (0, 0, 0) (m)')
plt.plot(index, objDist['objectDist'+str(0)], color = 'red', label = 'Object 0')
plt.plot(index, objDist['objectDist'+str(1)], color = 'green', label = 'Object 1')
plt.plot(index, objDist['objectDist'+str(2)], color = 'blue', label = 'Object 2')
plt.legend()
plt.show()
