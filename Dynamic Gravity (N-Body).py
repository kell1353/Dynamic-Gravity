from mayavi import mlab
import numpy as np
import math as m
import time

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
            a(u[i][1], u[i][4], u[i][5], u[i][6], xp, yp, zp)
            ax_sum += ax; ay_sum += ay; az_sum += az

'Compute the next x,y,z positions and velocities using the 4th-order Runge-Kutta technique'
def computeParameters(ID, x, y, z, vx, vy, vz):
    global x_n1; global y_n1; global z_n1; global vx_n1; global vy_n1; global vz_n1
    'Step Length: 1 day = 86,400 (in seconds (s))'
    del_t = 86400
    
    k_1x = vx 
    k_1y = vy 
    k_1z = vz
    a_sum(ID, x, y, z)
    k_1vx = ax_sum
    k_1vy = ay_sum
    k_1vz = az_sum
    
    k_2x = vx + ((del_t/2)*k_1vx)
    k_2y = vy + ((del_t/2)*k_1vy)
    k_2z = vz + ((del_t/2)*k_1vz)
    a_sum(ID, x + ((del_t/2)*k_1x), y + ((del_t/2)*k_1y), z + ((del_t/2)*k_1z))
    k_2vx = ax_sum
    k_2vy = ay_sum
    k_2vz = az_sum
    
    k_3x = vx + ((del_t/2)*k_2vx)
    k_3y = vy + ((del_t/2)*k_2vy)
    k_3z = vz + ((del_t/2)*k_2vz)
    a_sum(ID, x + ((del_t/2)*k_2x), y + ((del_t/2)*k_2y), z + ((del_t/2)*k_2z))
    k_3vx = ax_sum
    k_3vy = ay_sum
    k_3vz = az_sum

    k_4x = vx + ((del_t)*k_3vx)
    k_4y = vy + ((del_t)*k_3vy)
    k_4z = vz + ((del_t)*k_3vz)
    a_sum(ID, x + ((del_t)*k_3x), y + ((del_t)*k_3y), z + ((del_t)*k_3z))
    k_4vx = ax_sum
    k_4vy = ay_sum
    k_4vz = az_sum

    x_n1 = x + ((del_t/6)*(k_1x + (2*k_2x) + (2*k_3x) + k_4x))
    y_n1 = y + ((del_t/6)*(k_1y + (2*k_2y) + (2*k_3y) + k_4y))
    z_n1 = z + ((del_t/6)*(k_1z + (2*k_2z) + (2*k_3z) + k_4z))

    vx_n1 = vx + ((del_t/6)*(k_1vx + (2*k_2vx) + (2*k_3vx) + k_4vx))
    vy_n1 = vy + ((del_t/6)*(k_1vy + (2*k_2vy) + (2*k_3vy) + k_4vy))
    vz_n1 = vz + ((del_t/6)*(k_1vz + (2*k_2vz) + (2*k_3vz) + k_4vz))


#Add in an ID number, eq & polar radius, and a string for its name to grav objects
'Set Initial Positions (in meters (m)) and Set Initial Velocities (in m/s) and Mass (in kilograms (kg)) of the gravitational objects'
'gravObjects = [(ID, Mass, EqRadius, PolarRadius, Initial x_pos, Initial y_pos, Initial z_pos, Initial x_vel, Initial y_vel, Initial z_vel)]'   
#[0, 1.989*(10**30), 6955100000, 6955100000, -147095000000, 0, 0, 0, 0, 0]
u = [[0, 1.989*(10**30), 6955100000, 6955100000, 0, 0, 0, 0, 0, 0], \
         [1, 5.972*(10**24), 237100000, 237100000, 149600000000, 0, 0, 0, 30300, 0], \
         [2, 7.348*(10**22), 173710000, 173710000, 149215600000, 0, 0, 0, 31300, 0]]
#6371000 #1737100 #237100000 #384400000


'Set Initial Conditions'
'Gravitational Constant (in N*m^2/(kg^2))'
G = 6.67408*(10**(-11))
scaling_factor = 1550000000


'Draw Initial Objects'
'0'
draw_sphere(u[0][2], u[0][3])
sun = mlab.mesh(x, y, z)
'1'
draw_sphere(u[1][2], u[1][3])
earth = mlab.mesh(x, y, z)
'2'
draw_sphere(u[2][2], u[2][3])
moon = mlab.mesh(x, y, z)


days = 365
'Draw the Trajectory Lines for each objects updating for each objects position'
xtraj_0, ytraj_0, ztraj_0 = [], [], []
xtraj_1, ytraj_1, ztraj_1 = [], [], []
xtraj_2, ytraj_2, ztraj_2 = [], [], []
for i in range(0, days):
    for j in range(0, len(u)):
        computeParameters(u[j][0], u[j][4], u[j][5], u[j][6], u[j][7], u[j][8], u[j][9])
        if j == 0:
            xtraj_0.append(x_n1); ytraj_0.append(y_n1); ztraj_0.append(z_n1)
        elif j == 1:
            xtraj_1.append(x_n1); ytraj_1.append(y_n1); ztraj_1.append(z_n1)
        elif j == 2:
            xtraj_2.append(x_n1); ytraj_2.append(y_n1); ztraj_2.append(z_n1)
                
        u[j][4], u[j][5], u[j][6] = x_n1, y_n1, z_n1
        u[j][7], u[j][8], u[j][9] = vx_n1, vy_n1, vz_n1

mlab.plot3d(xtraj_0, ytraj_0, ztraj_0, color=(1,1,1), tube_radius = 100000000/2)
mlab.plot3d(xtraj_1, ytraj_1, ztraj_1, color=(1,1,1), tube_radius = 100000000/2)
mlab.plot3d(xtraj_2, ytraj_2, ztraj_2, color=(1,1,1), tube_radius = 75000000/2)


u = [[0, 1.989*(10**30), 6955100000, 6955100000, 0, 0, 0, 0, 0, 0], \
         [1, 5.972*(10**24), 237100000, 237100000, 149600000000, 0, 0, 0, 30300, 0], \
         [2, 7.348*(10**22), 173710000, 173710000, 149215600000, 0, 0, 0, 31300, 0]]

'Animate the figure'
@mlab.animate(delay = 150)
def anim():
    for i in range(0, days):
        time.sleep(.00000000000000001)
        for j in range(0, len(u)):
            computeParameters(u[j][0], u[j][4], u[j][5], u[j][6], u[j][7], u[j][8], u[j][9])
            if j == 0:
                sun.mlab_source.set(x = u[0][2] * np.sin(phi) * np.cos(theta) + x_n1, y = u[0][2] * np.sin(phi) * np.sin(theta) + y_n1, z = u[0][3] * np.cos(phi) + z_n1)
            elif j == 1:
                earth.mlab_source.set(x =  u[1][2] * np.sin(phi) * np.cos(theta) + x_n1, y = u[1][2] * np.sin(phi) * np.sin(theta) + y_n1, z = u[1][3] * np.cos(phi) + z_n1)
            elif j == 2:
                moon.mlab_source.set(x =  u[2][2] * np.sin(phi) * np.cos(theta) + x_n1, y = u[2][2] * np.sin(phi) * np.sin(theta) + y_n1, z = u[2][3] * np.cos(phi) + z_n1)
                
            u[j][4], u[j][5], u[j][6] = x_n1, y_n1, z_n1
            u[j][7], u[j][8], u[j][9] = vx_n1, vy_n1, vz_n1
        yield
    print(i)
anim()
mlab.show()
