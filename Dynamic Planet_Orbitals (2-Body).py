from mayavi import mlab
import numpy as np
import math as m
import time

fig = mlab.figure('Solar System', bgcolor = (0,0,0), size = (700,500))

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

def a_sum(xp, yp, zp):
    global ax_sum; global ay_sum; global az_sum  
    ax_sum, ay_sum, az_sum = 0, 0, 0
    for i in range(len(gravObjects)):
        a(gravObjects[i][0], gravObjects[i][1], gravObjects[i][2], gravObjects[i][3], xp, yp, zp)
        ax_sum += ax
        ay_sum += ay
        az_sum += az
        
'Compute the next x,y,z positions and velocities using the 4th-order Runge-Kutta technique'
def computeParameters(x, y, z, vx, vy, vz):
    global x_n1; global y_n1; global z_n1; global vx_n1; global vy_n1; global vz_n1
    'Step Length: 1 day = 86,400 (in seconds (s))'
    del_t = 86400
    
    k_1x = vx 
    k_1y = vy 
    k_1z = vz
    a_sum(x, y, z)
    k_1vx = ax_sum
    k_1vy = ay_sum
    k_1vz = az_sum
    
    k_2x = vx + ((del_t/2)*k_1vx)
    k_2y = vy + ((del_t/2)*k_1vy)
    k_2z = vz + ((del_t/2)*k_1vz)
    a_sum(x + ((del_t/2)*k_1x), y + ((del_t/2)*k_1y), z + ((del_t/2)*k_1z))
    k_2vx = ax_sum
    k_2vy = ay_sum
    k_2vz = az_sum
    
    k_3x = vx + ((del_t/2)*k_2vx)
    k_3y = vy + ((del_t/2)*k_2vy)
    k_3z = vz + ((del_t/2)*k_2vz)
    a_sum(x + ((del_t/2)*k_2x), y + ((del_t/2)*k_2y), z + ((del_t/2)*k_2z))
    k_3vx = ax_sum
    k_3vy = ay_sum
    k_3vz = az_sum

    k_4x = vx + ((del_t)*k_3vx)
    k_4y = vy + ((del_t)*k_3vy)
    k_4z = vz + ((del_t)*k_3vz)
    a_sum(x + ((del_t)*k_3x), y + ((del_t)*k_3y), z + ((del_t)*k_3z))
    k_4vx = ax_sum
    k_4vy = ay_sum
    k_4vz = az_sum

    x_n1 = x + ((del_t/6)*(k_1x + (2*k_2x) + (2*k_3x) + k_4x))
    y_n1 = y + ((del_t/6)*(k_1y + (2*k_2y) + (2*k_3y) + k_4y))
    z_n1 = z + ((del_t/6)*(k_1z + (2*k_2z) + (2*k_3z) + k_4z))
    
    vx_n1 = vx + ((del_t/6)*(k_1vx + (2*k_2vx) + (2*k_3vx) + k_4vx))
    vy_n1 = vy + ((del_t/6)*(k_1vy + (2*k_2vy) + (2*k_3vy) + k_4vy))
    vz_n1 = vz + ((del_t/6)*(k_1vz + (2*k_2vz) + (2*k_3vz) + k_4vz))
    

'gravObjects = [(Mass (in kilograms (kg)), Initial x pos, Initial y pos, Initial z pos), .....]'   
##gravObjects = [(1.989*(10**30), -147095000000, -147095000000, 0), \
##                              (1.989*(10**30), 147095000000, 147095000000, 0)]
gravObjects = [[1.989*(10**30), 0, 0, 0]]
'Set Initial Conditions'
'Gravitational Constant (in N*m^2/(kg^2))'
G = 6.67408*(10**(-11))
scaling_factor = 1550000000

'Stellar Constants'
'Stellar Radius (in m)'
R = 69595600*50
'Set Initial Positions (in meters (m))'
xc1, yc1, zc1 = gravObjects[0][1], gravObjects[0][2], gravObjects[0][3]
'Set Initial Velocities (in m/s)'
vx_c, vy_c, vz_c = 0, 0, 0
'Draw Stars'
draw_object(R, xc1, yc1, zc1)


'Planetary Constants'
'Set Object Masses (in kilograms (kg))'
m = 5.972*(10**24)
'Set Object Radii (in meters (m))'
me_eqRad, me_polarRad = 243970000, 243970000
v_eqRad, v_polarRad = 605180000, 605180000
e_eqRad, e_polarRad = 637810000, 635680000
ma_eqRad, ma_polarRad = 339620000, 337620000
j_eqRad, j_polarRad = 7149200000, 6685400000
s_eqRad, s_polarRad = 6026800000, 5436400000
u_eqRad, u_polarRad = 2555900000, 2497300000
n_eqRad, n_polarRad = 2476400000, 2434100000
'Set Initial Positions (in meters (m))'
me_xp0, me_yp0, me_zp0 = -57910000000, 0, 0   
v_xp0, v_yp0, v_zp0 = -108210000000, 0, 0
e_xp0, e_yp0, e_zp0 = -147095000000, 0, 0
ma_xp0, ma_yp0, ma_zp0 = -227920000000, 0, 0
j_xp0, j_yp0, j_zp0 = -778570000000, 0, 0
s_xp0, s_yp0, s_zp0 = -1433530000000, 0, 0
u_xp0, u_yp0, u_zp0 = -2872460000000, 0, 0
n_xp0, n_yp0, n_zp0 = -5906380000000, 0, 0
'Set Initial Velocities (in m/s)'
me_vxp0, me_vyp0, me_vzp0 = 0, -47340, 0           
v_vxp0, v_vyp0, v_vzp0 = 0, -35012, 0
e_vxp0, e_vyp0, e_vzp0 = 0, -30300, 0
ma_vxp0, ma_vyp0, ma_vzp0 = 0, -24072, 0
j_vxp0, j_vyp0, j_vzp0 = 0, -13040, 0
s_vxp0, s_vyp0, s_vzp0 = 0, -9611, 0
u_vxp0, u_vyp0, u_vzp0 = 0, -5431, 0
n_vxp0, n_vyp0, n_vzp0 = 0, -7006, 0 
'Draw Initial Object'
draw_sphere(me_eqRad, me_polarRad)
mercury = mlab.mesh(x, y, z)
draw_sphere(v_eqRad, v_polarRad)
venus = mlab.mesh(x, y, z)
draw_sphere(e_eqRad, e_polarRad)
earth = mlab.mesh(x, y, z)
draw_sphere(ma_eqRad, ma_polarRad)
mars = mlab.mesh(x, y, z)
draw_sphere(j_eqRad, j_polarRad)
jupiter = mlab.mesh(x, y, z)
draw_sphere(s_eqRad, s_polarRad)
saturn = mlab.mesh(x, y, z)
draw_sphere(u_eqRad, u_polarRad)
uranus = mlab.mesh(x, y, z)
draw_sphere(n_eqRad, n_polarRad)
neptune = mlab.mesh(x, y, z)
#velocity_text = mlab.text3d(-147095000000, -147095000000, 0 + (polarRad*2), str(vectorMag(vx_p0, vy_p0, vz_p0)), scale=(scaling_factor, scaling_factor, scaling_factor))


'Calculate and Plot the Trajectory Line'
def calcTraj(d, rad, x, y, z ,vx, vy, vz):
    x_traj , y_traj , z_traj = [], [], []
    for i in range(0, d):
        if i == 1:    
            computeParameters(x, y, z , vx, vy, vz)
            x_traj.append(x_n1),    y_traj.append(y_n1),    z_traj.append(z_n1)
        elif i > 1:
            computeParameters(x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1)
            x_traj.append(x_n1),    y_traj.append(y_n1),    z_traj.append(z_n1)
    mlab.plot3d(x_traj , y_traj , z_traj, color=(1,1,1), tube_radius = rad)



days = 124
calcTraj(days, 100000000/3, me_xp0, me_yp0, me_zp0, me_vxp0, me_vyp0, me_vzp0)
calcTraj(days, 100000000/2, v_xp0, v_yp0, v_zp0, v_vxp0, v_vyp0, v_vzp0)
calcTraj(days, 100000000/2, e_xp0, e_yp0, e_zp0, e_vxp0, e_vyp0, e_vzp0)
calcTraj(days, 100000000/2, ma_xp0, ma_yp0, ma_zp0, ma_vxp0, ma_vyp0, ma_vzp0)
calcTraj(days, 100000000, j_xp0, j_yp0, j_zp0, j_vxp0, j_vyp0, j_vzp0)
calcTraj(days, 100000000, s_xp0, s_yp0, s_zp0, s_vxp0, s_vyp0, s_vzp0)
calcTraj(days, 100000000, u_xp0, u_yp0, u_zp0, u_vxp0, u_vyp0, u_vzp0)
calcTraj(days, 100000000, n_xp0, n_yp0, n_zp0, n_vxp0, n_vyp0, n_vzp0)


'Animate the figure'
@mlab.animate(delay = 10)
def anim():
    for i in range(0, days):
        time.sleep(.000000000001)
        if i == 1:
            'Mercury'
            computeParameters(me_xp0, me_yp0, me_zp0, me_vxp0, me_vyp0, me_vzp0)
            me_x, me_y, me_z, me_vx, me_vy, me_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            mercury.mlab_source.set(x = me_eqRad * np.sin(phi) * np.cos(theta) + me_x, y = me_eqRad * np.sin(phi) * np.sin(theta) + me_y, z = me_polarRad * np.cos(phi) + me_z)
            'Venus'
            computeParameters(v_xp0, v_yp0, v_zp0, v_vxp0, v_vyp0, v_vzp0)
            v_x, v_y, v_z, v_vx, v_vy, v_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            venus.mlab_source.set(x = v_eqRad * np.sin(phi) * np.cos(theta) + v_x, y = v_eqRad * np.sin(phi) * np.sin(theta) + v_y, z = v_polarRad * np.cos(phi) + v_z)
            'Earth'
            computeParameters(e_xp0, e_yp0, e_zp0, e_vxp0, e_vyp0, e_vzp0)
            e_x, e_y, e_z, e_vx, e_vy, e_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            earth.mlab_source.set(x = e_eqRad * np.sin(phi) * np.cos(theta) + e_x, y = e_eqRad * np.sin(phi) * np.sin(theta) + e_y, z = e_polarRad * np.cos(phi) + e_z)
            'Mars'
            computeParameters(ma_xp0, ma_yp0, ma_zp0, ma_vxp0, ma_vyp0, ma_vzp0)
            ma_x, ma_y, ma_z, ma_vx, ma_vy, ma_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            mars.mlab_source.set(x = ma_eqRad * np.sin(phi) * np.cos(theta) + ma_x, y = e_eqRad * np.sin(phi) * np.sin(theta) + ma_y, z = e_polarRad * np.cos(phi) + ma_z)
            'Jupiter'
            computeParameters(j_xp0, j_yp0, j_zp0, j_vxp0, j_vyp0, j_vzp0)
            j_x, j_y, j_z, j_vx, j_vy, j_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            jupiter.mlab_source.set(x = j_eqRad * np.sin(phi) * np.cos(theta) + j_x, y = j_eqRad * np.sin(phi) * np.sin(theta) + j_y, z = j_polarRad * np.cos(phi) + j_z)
            'Saturn'
            computeParameters(s_xp0, s_yp0, s_zp0, s_vxp0, s_vyp0, s_vzp0)
            s_x, s_y, s_z, s_vx, s_vy, s_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            saturn.mlab_source.set(x = s_eqRad * np.sin(phi) * np.cos(theta) + s_x, y = s_eqRad * np.sin(phi) * np.sin(theta) + s_y, z = s_polarRad * np.cos(phi) + s_z)
            'Uranus'
            computeParameters(u_xp0, u_yp0, u_zp0, u_vxp0, u_vyp0, u_vzp0)
            u_x, u_y, u_z, u_vx, u_vy, u_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            uranus.mlab_source.set(x = u_eqRad * np.sin(phi) * np.cos(theta) + u_x, y = u_eqRad * np.sin(phi) * np.sin(theta) + u_y, z = u_polarRad * np.cos(phi) + u_z)
            'Neptune'
            computeParameters(n_xp0, n_yp0, n_zp0, n_vxp0, n_vyp0, n_vzp0)
            n_x, n_y, n_z, n_vx, n_vy, n_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            neptune.mlab_source.set(x = n_eqRad * np.sin(phi) * np.cos(theta) + n_x, y = n_eqRad * np.sin(phi) * np.sin(theta) + n_y, z = n_polarRad * np.cos(phi) + n_z)
        elif i > 1:
            'Mercury'
            computeParameters(me_x, me_y, me_z, me_vx, me_vy, me_vz)
            me_x, me_y, me_z, me_vx, me_vy, me_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            mercury.mlab_source.set(x = me_eqRad * np.sin(phi) * np.cos(theta) + me_x, y = me_eqRad * np.sin(phi) * np.sin(theta) + me_y, z = me_polarRad * np.cos(phi) + me_z)
            'Venus'
            computeParameters(v_x, v_y, v_z, v_vx, v_vy, v_vz)
            v_x, v_y, v_z, v_vx, v_vy, v_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            venus.mlab_source.set(x = v_eqRad * np.sin(phi) * np.cos(theta) + v_x, y = v_eqRad * np.sin(phi) * np.sin(theta) + v_y, z = v_polarRad * np.cos(phi) + v_z)
            'Earth'
            computeParameters(e_x, e_y, e_z, e_vx, e_vy, e_vz)
            e_x, e_y, e_z, e_vx, e_vy, e_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            earth.mlab_source.set(x = e_eqRad * np.sin(phi) * np.cos(theta) + e_x, y = e_eqRad * np.sin(phi) * np.sin(theta) + e_y, z = e_polarRad * np.cos(phi) + e_z)
            'Mars'
            computeParameters(ma_x, ma_y, ma_z, ma_vx, ma_vy, ma_vz)
            ma_x, ma_y, ma_z, ma_vx, ma_vy, ma_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            mars.mlab_source.set(x = ma_eqRad * np.sin(phi) * np.cos(theta) + ma_x, y = e_eqRad * np.sin(phi) * np.sin(theta) + ma_y, z = e_polarRad * np.cos(phi) + ma_z)
            'Jupiter'
            computeParameters(j_x, j_y, j_z, j_vx, j_vy, j_vz)
            j_x, j_y, j_z, j_vx, j_vy, j_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            jupiter.mlab_source.set(x = j_eqRad * np.sin(phi) * np.cos(theta) + j_x, y = j_eqRad * np.sin(phi) * np.sin(theta) + j_y, z = j_polarRad * np.cos(phi) + j_z)
            'Saturn'
            computeParameters(s_x, s_y, s_z, s_vx, s_vy, s_vz)
            s_x, s_y, s_z, s_vx, s_vy, s_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            saturn.mlab_source.set(x = s_eqRad * np.sin(phi) * np.cos(theta) + s_x, y = s_eqRad * np.sin(phi) * np.sin(theta) + s_y, z = s_polarRad * np.cos(phi) + s_z)
            'Uranus'
            computeParameters(u_x, u_y, u_z, u_vx, u_vy, u_vz)
            u_x, u_y, u_z, u_vx, u_vy, u_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            uranus.mlab_source.set(x = u_eqRad * np.sin(phi) * np.cos(theta) + u_x, y = u_eqRad * np.sin(phi) * np.sin(theta) + u_y, z = u_polarRad * np.cos(phi) + u_z)
            'Neptune'
            computeParameters(n_x, n_y, n_z, n_vx, n_vy, n_vz)
            n_x, n_y, n_z, n_vx, n_vy, n_vz = x_n1, y_n1, z_n1, vx_n1, vy_n1, vz_n1
            neptune.mlab_source.set(x = n_eqRad * np.sin(phi) * np.cos(theta) + n_x, y = n_eqRad * np.sin(phi) * np.sin(theta) + n_y, z = n_polarRad * np.cos(phi) + n_z)
        yield

anim()
mlab.show()
