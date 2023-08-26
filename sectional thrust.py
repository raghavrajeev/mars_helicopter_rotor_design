# -*- coding: utf-8 -*-

import math
import numpy as np
import matplotlib.pyplot as plt

#  Density, Radius, root radius, number of blades, chord and mass respectively
rho=0.02
R=0.5
Rrc=0.0525
b=3
c=0.0508
m=2
g=3.72

# RPM, angular velocity and vertical velocity
rpm=2000
omega = 2*math.pi*rpm/60

# Lift coefficient terms:
a=8

# # Drag coefficient terms:
# cd0=0.0113
# cd1=1.25

# Coefficient of drag (fuselage) and wetted area
Cd_f = 0.4
S=0.9

# Forward velocity
V_inf=30 #100 km/hr = 27.7 m/s
v=6

# alpha_tpp calculation
drag = 0.5 * rho * (Cd_f) * V_inf**2 * S
weight=m*g
alpha_tpp = math.atan(drag/weight)

# Thrust, thrust coefficient and advance ratio calculation
T = weight/math.cos(alpha_tpp)
C_T = T/(rho*math.pi*(R**2)*(omega*R)**2)
mu=(V_inf*math.cos(alpha_tpp))/(omega*R)


# Inflow model
lamb_ig = C_T/(mu*2)
lamb_g =lamb_ig+(V_inf*math.sin(alpha_tpp)/(omega*R))

def inflow(r, psi): # non-uniform main rotor inflow
    lamb_i = lamb_ig*(1+(((4*mu)/3)/(mu+1.2*(lamb_g)))*(r/R)*math.cos(psi))
    lamb_i = lamb_i + ((V_inf*math.sin(alpha_tpp))/(omega*R))
    return (lamb_i)

# Angles (Pitch and flap)
# collective pitch, lateral and longitudinal cyclic
theta_0 = 27
theta_tw = 12
theta_1c = 5.61
theta_1s = -12.43
# flap angle, longitudinal and lateral tilt
beta_0 = 24

# converting to radians
theta_0 = theta_0*math.pi/180
theta_tw = theta_tw*math.pi/180
theta_1s = theta_1s*math.pi/180
theta_1c = theta_1c*math.pi/180
beta_0 = beta_0*math.pi/180

def pitch(r, psi):
    return (theta_0 + theta_tw*(r/R) + theta_1c*math.cos(psi) + theta_1s*math.sin(psi))
# Function to calculate Pitch values. 

def aoa(r, psi):
    U_P = inflow(r,psi)*omega*R+V_inf*math.cos(psi)*math.sin(beta_0)+6
    U_T = omega*r + V_inf*math.sin(psi)*math.cos(alpha_tpp)
    return (pitch(r, psi) - math.atan(U_P/U_T))

rpts=[]
for i in range(2000):
    rpts.append(Rrc + i*(R-Rrc)/2000)
#Dividing rotor blade into 2000 segmets

psipts=[]
for i in range(360):
    psipts.append(i*math.pi/180)
# Dividing Azimuthal angles into 360 components

#Sectional lift
def lift(r,psi):
    t1=theta_0 * omega**2 * r**2 + 2* V_inf**2 * omega * r * theta_0 *math.sin(psi) + theta_0*V_inf**2 * (math.sin(psi))**2
    + theta_tw * omega**2 * r**3 / R + 2*V_inf * omega *theta_tw * r**2 * math.sin(psi) /R +theta_tw*(r/R)*(math.sin(psi))**2
    t2 = omega**2 * r**2 * theta_1c * math.cos(psi) + 2 * omega * r * V_inf * theta_1c * math.sin(psi) * math.cos(psi) + V_inf**2 * theta_1c * (math.sin(psi))**2 * math.cos(psi)
    + omega**2 + r**2 + theta_1s * math.sin(psi) + V_inf**2 + theta_1s * (math.sin(psi))**3 + 2 * omega * r * V_inf * theta_1s *(math.sin(psi))**2
    t3 = - omega**2 * R * inflow(r, psi) * r - V_inf * omega * R * inflow(r, psi) * math.sin(psi) 
    -omega * r * V_inf * beta_0 * math.cos(psi) - V_inf**2 * beta_0 * math.cos(psi) * math.sin(psi)
    L=t1+t2+t3
    return L
    
# Calculating thrust from sectional lift
def thrust(r,psi):
    T=b*(1/(2*math.pi))*lift(r,psi)
    return T

# Assigning coordinates for contour plot
def Rotate_for_plot(x2, y2, theta):
    X = (x2 * np.cos(theta) - y2 * np.sin(theta))
    Y = (x2 * np.sin(theta) + y2 * np.cos(theta))
    return X,Y

radius=np.linspace(Rrc,R,2000)

x_values= np.zeros((361,2000))  # N is number of azimuthal divisions , IM  is number of radial divisions
y_values= np.zeros((361,2000))  
z_values= np.zeros((361,2000))
for j in range(361):
    for i in range(2000):
        x_values[j][i],y_values[j][i] = Rotate_for_plot(0,radius[i]/R, -j*(math.pi/180)) # dpsi is (360/(N-1))*(pi/180)  ,in radians , R_tip is tip radius

fig, ax = plt.subplots(1, 1, dpi=120, figsize=(5,5))
# 'Thrust_plot_array' is a (N,IM) size array, rows represent azimuthal location and columns represent radial location
#  You should fill 'Thrust_plot_array' inside loop (according to your way of making the code)
Thrust_plot_array=np.zeros((360,2000))
for j in range(360):
    for i in range(2000):
        Thrust_plot_array[j][i]=thrust(radius[i], j*(math.pi/180))
for j in range(361):
    if j==360:
        z_values[j,:] = Thrust_plot_array [0,:]   
    else:
        z_values[j,:] = Thrust_plot_array [j,:]

ax.contourf(x_values,y_values, z_values)
ax.set_title('Thrust Contour')
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.show()